function factor_correlation_analysis(results_dir)
%FACTOR_CORRELATION_ANALYSIS Which per-epoch factors track the WLS horizontal error?
%
%   factor_correlation_analysis('results/samsung3rd_factors')
%
% Reads factor_epochs.csv (one row per solved epoch, written by the Python
% exporter: python -m wls_outlier_lab.cli ... --factor-export) and runs the
% correlation study in MATLAB:
%   * speed-interval segmentation (static / low / mid / high from truth speed)
%   * Pearson + Spearman correlation of every factor vs horizontal error,
%     overall and per speed bin
%   * partial Spearman controlling geometry (HDOP + n_used)
%   * moving-block bootstrap 95% CI on Spearman (30 s blocks — the error series
%     is autocorrelated, naive p-values would overstate significance)
%   * factor ranking, figures, CSV/TXT outputs
%
% Outputs (into results_dir):
%   factor_correlations.csv   long format: factor x group correlations + CI
%   factor_bin_stats.csv      horizontal-error stats per speed bin
%   factor_summary.txt        ranked findings (text)
%   matlab_factor_heatmap.png, matlab_factor_rank.png,
%   matlab_factor_top_scatter.png, matlab_factor_timeseries.png
%
% No Statistics Toolbox required: rank / correlation / quantile / bootstrap are
% implemented locally.

    if nargin < 1, results_dir = 'results/samsung3rd_factors'; end
    fe_file = fullfile(results_dir, 'factor_epochs.csv');
    assert(isfile(fe_file), 'factor_epochs.csv not found in %s (run the Python exporter first)', results_dir);
    T = readtable(fe_file, 'VariableNamingRule', 'preserve');
    fprintf('[factor] %d epochs x %d columns\n', height(T), width(T));

    rng(42);  % reproducible bootstrap

    % --- palette (validated defaults; see dataviz notes in the repo report) --
    C.ink     = [11 11 11]/255;
    C.muted   = [137 135 129]/255;
    C.grid    = [225 224 217]/255;
    C.blue    = [42 120 214]/255;     % categorical slot 1
    C.orange  = [235 104 52]/255;     % categorical slot 2
    C.aqua    = [27 175 122]/255;     % categorical slot 3
    C.red     = [227 73 72]/255;      % categorical slot 8 (kept for h_err, matches repo figures)
    C.seq     = [134 182 239; 57 135 229; 28 92 171; 13 54 107]/255;  % ordinal blues 250/400/550/700
    C.divlo   = [28 92 171]/255;      % diverging blue pole
    C.divmid  = [240 239 236]/255;    % neutral midpoint
    C.divhi   = [208 59 59]/255;      % diverging red pole

    t   = T.t_sec - min(T.t_sec);
    h   = T.h_err_m;
    spd = T.speed_mps;

    % --- derived factors --------------------------------------------------
    % clk_m drifts ~-242 m/s, so raw clk_m vs error is only a time-trend proxy;
    % the detrended bias and the drift *deviations* are the meaningful forms.
    fin = isfinite(T.clk_m) & isfinite(t);
    p = polyfit(t(fin), T.clk_m(fin), 1);
    T.clk_detrend_m = T.clk_m - polyval(p, t);
    T.clk_drift_dev_mps  = abs(T.clk_drift_mps  - median(T.clk_drift_mps(isfinite(T.clk_drift_mps))));
    T.drift_dop_dev_mps  = abs(T.drift_dop_mps  - median(T.drift_dop_mps(isfinite(T.drift_dop_mps))));
    T.abs_accel_mps2 = abs(T.accel_mps2);

    excl = {'t_sec', 'dt_s', 'h_err_m', 'v_err_m'};
    names = T.Properties.VariableNames;
    factors = names(~ismember(names, excl));

    % --- speed bins --------------------------------------------------------
    edges  = [0 0.5 3 8 inf];
    glabel = {'all', 'static(<0.5)', 'low(0.5-3)', 'mid(3-8)', 'high(>8)'};
    gmask  = cell(1, 5);
    gmask{1} = true(height(T), 1);
    for b = 1:4
        gmask{b+1} = isfinite(spd) & spd >= edges(b) & spd < edges(b+1);
    end

    % --- bin-level error statistics ----------------------------------------
    fid = fopen(fullfile(results_dir, 'factor_bin_stats.csv'), 'w');
    fprintf(fid, 'bin,n,h_median_m,h_p95_m,h_rmse_m,speed_mean_mps\n');
    binstat = zeros(5, 5);
    for g = 1:5
        hh = h(gmask{g} & isfinite(h));
        ss = spd(gmask{g} & isfinite(spd));
        binstat(g, :) = [numel(hh), median(hh), local_quantile(hh, 0.95), ...
                         sqrt(mean(hh.^2)), mean(ss)];
        fprintf(fid, '%s,%d,%.3f,%.3f,%.3f,%.3f\n', glabel{g}, binstat(g,1), ...
                binstat(g,2), binstat(g,3), binstat(g,4), binstat(g,5));
    end
    fclose(fid);

    % --- correlations: factor x group --------------------------------------
    nf = numel(factors);
    MIN_N = 60;         % below this a bin correlation is too noisy to rank
    B_BOOT = 1000; L_BLOCK = 30;
    [R_pear, R_spear, R_part, R_part2, CI_lo, CI_hi, N_used] = deal(nan(nf, 5));
    ctrl_names = {'hdop', 'n_used', 'speed_mps'};
    ctrl = [T.hdop, T.n_used, T.speed_mps];
    for f = 1:nf
        x = T.(factors{f});
        if ~isnumeric(x), continue; end
        for g = 1:5
            m = gmask{g} & isfinite(x) & isfinite(h);
            n = nnz(m);
            N_used(f, g) = n;
            if n < MIN_N || numel(unique(x(m))) < 3, continue; end
            xs = x(m); hs = h(m);
            R_pear(f, g)  = local_pearson(xs, hs);
            R_spear(f, g) = local_spearman(xs, hs);
            keepc = ~strcmp(ctrl_names, factors{f});
            cc = ctrl(m, keepc & [true true false]);
            R_part(f, g)  = local_partial_spearman(hs, xs, cc);          % geometry only
            cc2 = ctrl(m, keepc);
            R_part2(f, g) = local_partial_spearman(hs, xs, cc2);         % geometry + speed
            % adaptive block length: keep >=6 blocks so small bins still get a CI
            Lg = min(L_BLOCK, max(8, floor(n / 6)));
            [CI_lo(f, g), CI_hi(f, g)] = local_block_boot(xs, hs, Lg, B_BOOT);
        end
    end

    % --- long-format CSV ----------------------------------------------------
    fid = fopen(fullfile(results_dir, 'factor_correlations.csv'), 'w');
    fprintf(fid, ['factor,bin,n,pearson_r,spearman_r,partial_geo_r,' ...
                  'partial_geo_speed_r,boot_lo,boot_hi,significant\n']);
    for f = 1:nf
        for g = 1:5
            sig = isfinite(CI_lo(f,g)) && (CI_lo(f,g) > 0 || CI_hi(f,g) < 0);
            fprintf(fid, '%s,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d\n', factors{f}, glabel{g}, ...
                    N_used(f,g), R_pear(f,g), R_spear(f,g), R_part(f,g), R_part2(f,g), ...
                    CI_lo(f,g), CI_hi(f,g), sig);
        end
    end
    fclose(fid);

    % --- ranking ------------------------------------------------------------
    [~, order] = sort(abs(R_spear(:, 1)), 'descend', 'MissingPlacement', 'last');
    sig1 = isfinite(CI_lo(:,1)) & (CI_lo(:,1) > 0 | CI_hi(:,1) < 0);

    fid = fopen(fullfile(results_dir, 'factor_summary.txt'), 'w');
    fprintf(fid, '=== Factor vs horizontal error: Spearman ranking (all epochs) ===\n');
    fprintf(fid, '%-22s %8s %8s %8s %8s %8s %19s %4s\n', 'factor', 'n', 'spear', 'pears', ...
            'p.geo', 'p.g+v', 'boot95CI', 'sig');
    for k = 1:nf
        f = order(k);
        if ~isfinite(R_spear(f, 1)), continue; end
        fprintf(fid, '%-22s %8d %8.3f %8.3f %8.3f %8.3f   [%6.3f, %6.3f] %4s\n', ...
                factors{f}, N_used(f,1), R_spear(f,1), R_pear(f,1), R_part(f,1), R_part2(f,1), ...
                CI_lo(f,1), CI_hi(f,1), tern(sig1(f), '*', ''));
    end
    fprintf(fid, ['\n(p.geo = partial r controlling HDOP+n_used; p.g+v also controls speed;\n' ...
                  ' sig ''*'' = 95%% moving-block bootstrap CI excludes 0; block <= %d s, B=%d)\n'], ...
            L_BLOCK, B_BOOT);
    fprintf(fid, '\n=== Per speed bin: top-5 |spearman| ===\n');
    for g = 2:5
        fprintf(fid, '\n-- %s  (n=%d, h_med=%.2f m, h_p95=%.2f m)\n', ...
                glabel{g}, binstat(g,1), binstat(g,2), binstat(g,3));
        if ~any(isfinite(R_spear(:, g)))
            fprintf(fid, '   (suppressed: n < %d epochs)\n', MIN_N);
            continue;
        end
        [~, og] = sort(abs(R_spear(:, g)), 'descend', 'MissingPlacement', 'last');
        for k = 1:min(5, nf)
            f = og(k);
            if ~isfinite(R_spear(f, g)), break; end
            sg = isfinite(CI_lo(f,g)) && (CI_lo(f,g) > 0 || CI_hi(f,g) < 0);
            fprintf(fid, '   %-22s r=%6.3f  partial=%6.3f  [%6.3f, %6.3f]%s\n', ...
                    factors{f}, R_spear(f,g), R_part(f,g), CI_lo(f,g), CI_hi(f,g), tern(sg, ' *', ''));
        end
    end
    fclose(fid);
    type(fullfile(results_dir, 'factor_summary.txt'));

    % =====================================================================
    % Figure 1: correlation heatmap (factors x groups), diverging colormap
    % =====================================================================
    keep = find(any(isfinite(R_spear), 2));
    [~, ord] = sort(abs(R_spear(keep, 1)), 'descend', 'MissingPlacement', 'last');
    rows = keep(ord);
    Rh = R_spear(rows, :);
    f1 = figure('Color', 'w', 'Name', 'Factor correlation heatmap', ...
                'Position', [40 40 860 30*numel(rows)+140], 'Visible', 'off');
    im = imagesc(Rh, [-0.7 0.7]);
    set(im, 'AlphaData', isfinite(Rh));            % NaN cells -> axes background
    set(gca, 'Color', [0.97 0.97 0.96]);
    colormap(f1, local_divmap(C.divlo, C.divmid, C.divhi));
    cb = colorbar; cb.Label.String = 'Spearman r vs horizontal error';
    set(gca, 'XTick', 1:5, 'XTickLabel', glabel, 'YTick', 1:numel(rows), ...
        'YTickLabel', strrep(factors(rows), '_', '\_'), 'TickLength', [0 0], ...
        'FontSize', 9, 'XColor', C.ink, 'YColor', C.ink);
    for i = 1:numel(rows)
        for g = 1:5
            v = Rh(i, g);
            if ~isfinite(v)
                text(g, i, '-', 'HorizontalAlignment', 'center', 'Color', C.muted, 'FontSize', 8);
                continue;
            end
            s = isfinite(CI_lo(rows(i),g)) && (CI_lo(rows(i),g) > 0 || CI_hi(rows(i),g) < 0);
            ink = tern(abs(v) > 0.45, [1 1 1], C.ink);
            text(g, i, sprintf('%.2f%s', v, tern(s, '*', '')), 'HorizontalAlignment', 'center', ...
                 'Color', ink, 'FontSize', 8, 'FontWeight', tern(s, 'bold', 'normal'));
        end
    end
    title('Factor vs horizontal error - Spearman r by speed interval (* = bootstrap 95% CI excludes 0)', ...
          'FontSize', 10, 'Color', C.ink);
    exportgraphics(f1, fullfile(results_dir, 'matlab_factor_heatmap.png'), 'Resolution', 130);

    % =====================================================================
    % Figure 2: signed ranking with bootstrap CI + error by speed bin
    % =====================================================================
    f2 = figure('Color', 'w', 'Name', 'Factor ranking', 'Position', [60 60 1240 560], 'Visible', 'off');
    subplot(1, 2, 1); hold on;
    ntop = min(14, numel(rows));
    rf = rows(1:ntop);
    vals = R_spear(rf, 1);
    for k = 1:ntop
        col = tern(vals(k) >= 0, C.divhi, C.divlo);
        barh(ntop - k + 1, vals(k), 0.62, 'FaceColor', col, 'EdgeColor', 'none');
        plot([CI_lo(rf(k),1) CI_hi(rf(k),1)], [ntop-k+1 ntop-k+1], '-', 'Color', C.ink, 'LineWidth', 1.1);
    end
    xline(0, '-', 'Color', C.muted);
    set(gca, 'YTick', 1:ntop, 'YTickLabel', strrep(flip(factors(rf)), '_', '\_'), ...
        'FontSize', 9, 'XColor', C.ink, 'YColor', C.ink, 'Box', 'off');
    grid on; set(gca, 'GridColor', C.grid, 'GridAlpha', 1);
    xlabel('Spearman r (all epochs), whisker = block-bootstrap 95% CI');
    title('Top factors by |r|  (red = error grows with factor)', 'FontSize', 10);

    subplot(1, 2, 2); hold on;
    xb = 1:4;
    b1 = bar(xb - 0.18, binstat(2:5, 2), 0.32, 'FaceColor', C.blue,   'EdgeColor', 'none');
    b2 = bar(xb + 0.18, binstat(2:5, 3), 0.32, 'FaceColor', C.orange, 'EdgeColor', 'none');
    for g = 1:4
        text(xb(g)-0.18, binstat(g+1,2)+0.08, sprintf('%.2f', binstat(g+1,2)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8.5, 'Color', C.ink);
        text(xb(g)+0.18, binstat(g+1,3)+0.08, sprintf('%.2f', binstat(g+1,3)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8.5, 'Color', C.ink);
        text(xb(g), -0.45, sprintf('n=%d', binstat(g+1,1)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', C.muted);
    end
    set(gca, 'XTick', xb, 'XTickLabel', glabel(2:5), 'FontSize', 9, ...
        'XColor', C.ink, 'YColor', C.ink, 'Box', 'off');
    grid on; set(gca, 'GridColor', C.grid, 'GridAlpha', 1);
    ylabel('horizontal error [m]');
    legend([b1 b2], {'median', 'p95'}, 'Location', 'northwest', 'Box', 'off');
    title('Horizontal error by speed interval', 'FontSize', 10);
    exportgraphics(f2, fullfile(results_dir, 'matlab_factor_rank.png'), 'Resolution', 130);

    % =====================================================================
    % Figure 3: top-6 factor scatters, points shaded by speed bin (ordinal)
    % =====================================================================
    % clk_m stays in the heatmap but is excluded here: it is monotone in time,
    % so its scatter only restates the session's error trend.
    cand = rows(~strcmp(factors(rows), 'clk_m'));
    ntop = min(6, numel(cand));
    f3 = figure('Color', 'w', 'Name', 'Top factor scatter', 'Position', [60 60 1380 760], 'Visible', 'off');
    for k = 1:ntop
        subplot(2, 3, k); hold on;
        fx = T.(factors{cand(k)});
        for g = 2:5
            m = gmask{g} & isfinite(fx) & isfinite(h);
            scatter(fx(m), h(m), 7, C.seq(g-1, :), 'filled', 'MarkerFaceAlpha', 0.45);
        end
        m = isfinite(fx) & isfinite(h);
        [xq, yq] = local_binned_median(fx(m), h(m), 10);
        plot(xq, yq, '-', 'Color', C.orange, 'LineWidth', 2);
        grid on; set(gca, 'GridColor', C.grid, 'GridAlpha', 1, 'FontSize', 8.5, ...
            'XColor', C.ink, 'YColor', C.ink, 'Box', 'off');
        xlabel(strrep(factors{cand(k)}, '_', '\_'));
        ylabel('horizontal error [m]');
        title(sprintf('r_s = %.2f   partial = %.2f', R_spear(cand(k),1), R_part(cand(k),1)), 'FontSize', 9.5);
        if k == 1
            legend([arrayfun(@(g) scatter(nan, nan, 12, C.seq(g,:), 'filled'), 1:4), ...
                    plot(nan, nan, '-', 'Color', C.orange, 'LineWidth', 2)], ...
                   [glabel(2:5), {'decile median'}], 'Location', 'best', 'Box', 'off', 'FontSize', 7.5);
        end
    end
    sgtitle('Top factors vs horizontal error (shade = speed interval)', 'FontSize', 11, 'Color', C.ink);
    exportgraphics(f3, fullfile(results_dir, 'matlab_factor_top_scatter.png'), 'Resolution', 130);

    % =====================================================================
    % Figure 4: time series - error, speed, two leading factors
    % =====================================================================
    lead = cand(1:min(2, numel(cand)));
    f4 = figure('Color', 'w', 'Name', 'Factor time series', 'Position', [60 60 1240 820], 'Visible', 'off');
    panels = [{h, 'horizontal error [m]', C.red}; {spd, 'truth speed [m/s]', C.blue}];
    for k = 1:numel(lead)
        panels(end+1, :) = {T.(factors{lead(k)}), strrep(factors{lead(k)}, '_', '\_'), ...
                            tern(k == 1, C.aqua, C.orange)}; %#ok<AGROW>
    end
    np = size(panels, 1);
    ax = gobjects(np, 1);
    for k = 1:np
        ax(k) = subplot(np, 1, k);
        plot(t, panels{k, 1}, '-', 'Color', panels{k, 3}, 'LineWidth', 1.0);
        ylabel(panels{k, 2}, 'FontSize', 9);
        grid on; set(gca, 'GridColor', C.grid, 'GridAlpha', 1, 'FontSize', 8.5, ...
            'XColor', C.ink, 'YColor', C.ink, 'Box', 'off');
        if k < np, set(gca, 'XTickLabel', []); end
    end
    linkaxes(ax, 'x'); xlim(ax(1), [min(t) max(t)]);
    xlabel(ax(end), 'time since start [s]');
    sgtitle('Horizontal error, speed and the two leading factors', 'FontSize', 11, 'Color', C.ink);
    exportgraphics(f4, fullfile(results_dir, 'matlab_factor_timeseries.png'), 'Resolution', 130);

    close([f1 f2 f3 f4]);
    fprintf('[factor] wrote factor_correlations.csv, factor_bin_stats.csv, factor_summary.txt + 4 PNGs to %s\n', results_dir);
end

% ========================= local statistics =============================

function r = local_rank(x)
%LOCAL_RANK average ranks with ties (input must be finite).
    n = numel(x);
    [sx, order] = sort(x(:));
    rr = zeros(n, 1);
    i = 1;
    while i <= n
        j = i;
        while j < n && sx(j+1) == sx(i), j = j + 1; end
        rr(i:j) = (i + j) / 2;
        i = j + 1;
    end
    r = zeros(n, 1);
    r(order) = rr;
end

function r = local_pearson(a, b)
    a = a(:); b = b(:);
    if numel(a) < 8 || std(a) == 0 || std(b) == 0, r = nan; return; end
    a = a - mean(a); b = b - mean(b);
    r = (a' * b) / sqrt((a' * a) * (b' * b));
end

function r = local_spearman(a, b)
    r = local_pearson(local_rank(a), local_rank(b));
end

function r = local_partial_spearman(y, x, ctrl)
%LOCAL_PARTIAL_SPEARMAN rank-transform, regress controls out of both, correlate.
    m = isfinite(y) & isfinite(x) & all(isfinite(ctrl), 2);
    if nnz(m) < 20, r = nan; return; end
    y = local_rank(y(m)); x = local_rank(x(m));
    Cr = ctrl(m, :);
    for c = 1:size(Cr, 2), Cr(:, c) = local_rank(Cr(:, c)); end
    A = [ones(nnz(m), 1), Cr];
    ry = y - A * (A \ y);
    rx = x - A * (A \ x);
    r = local_pearson(rx, ry);
end

function [lo, hi] = local_block_boot(x, y, L, B)
%LOCAL_BLOCK_BOOT moving-block bootstrap 95% CI of Spearman r.
    n = numel(x);
    lo = nan; hi = nan;
    if n < 3 * L, return; end
    nb = ceil(n / L);
    rs = nan(B, 1);
    for b = 1:B
        starts = randi(n - L + 1, nb, 1);
        idx = reshape(starts' + (0:L-1)', [], 1);
        idx = idx(1:n);
        rs(b) = local_spearman(x(idx), y(idx));
    end
    rs = sort(rs(isfinite(rs)));
    if numel(rs) < B / 2, return; end
    lo = rs(max(1, round(0.025 * numel(rs))));
    hi = rs(round(0.975 * numel(rs)));
end

function q = local_quantile(v, p)
    v = sort(v(isfinite(v)));
    if isempty(v), q = nan; return; end
    q = v(max(1, min(numel(v), round(p * numel(v)))));
end

function [xc, yc] = local_binned_median(x, y, nbin)
%LOCAL_BINNED_MEDIAN robust trend: median y within x-deciles.
    e = local_quantiles(x, nbin);
    xc = nan(nbin, 1); yc = nan(nbin, 1);
    for i = 1:nbin
        m = x >= e(i) & x <= e(i + 1);
        if nnz(m) < 5, continue; end
        xc(i) = median(x(m)); yc(i) = median(y(m));
    end
    m = isfinite(xc); xc = xc(m); yc = yc(m);
end

function e = local_quantiles(x, nbin)
    s = sort(x(:));
    e = s(max(1, round(linspace(1, numel(s), nbin + 1))));
    e(1) = e(1) - eps(abs(e(1)) + 1); e(end) = e(end) + eps(abs(e(end)) + 1);
end

function map = local_divmap(lo, mid, hi)
    n = 32;
    map = [interp1([0 1], [lo; mid], linspace(0, 1, n)); ...
           interp1([0 1], [mid; hi], linspace(0, 1, n))];
end

function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end
