function plot_wls_results(output_dir)
%PLOT_WLS_RESULTS  wls_outlier_lab 결과 CSV를 읽어 오차/비교 그림을 그린다.
%
%   plot_wls_results(OUTPUT_DIR)
%
% OUTPUT_DIR 는 파이썬 파이프라인(python -m wls_outlier_lab.cli ... --output-dir)
% 이 만든 폴더로, 다음 CSV들을 포함한다:
%   detector_comparison.csv      검출기별 오차 요약
%   per_epoch_<detector>.csv     에폭별 ENU 오차 (best/baseline)
%   constellation_ablation.csv   위성군 ablation (선택)
%
% 계산(WLS/이상치/메트릭)은 파이썬이 하고, 시각화만 여기서 담당한다.
% 그림은 화면에 뜨고 OUTPUT_DIR 에 matlab_*.png 로도 저장된다.

    if nargin < 1 || isempty(output_dir)
        error('plot_wls_results:usage', 'output_dir 를 지정하세요.');
    end
    assert(isfolder(output_dir), 'output_dir 가 없습니다: %s', output_dir);

    RED  = [0.753 0.224 0.169];
    BLUE = [0.173 0.435 0.733];

    comp = readtable(fullfile(output_dir, 'detector_comparison.csv'), ...
                     'VariableNamingRule', 'preserve');
    names = string(comp.detector);
    hrmse_all = comp.horizontal_rmse_m;
    [~, ibest] = min(hrmse_all);
    best = char(names(ibest));

    % ---- 에폭별 오차 (best 검출기) ---------------------------------------
    pe_file = fullfile(output_dir, sprintf('per_epoch_%s.csv', best));
    assert(isfile(pe_file), 'per-epoch 파일이 없습니다: %s', pe_file);
    pe = readtable(pe_file, 'VariableNamingRule', 'preserve');
    E = pe.east_m;  N = pe.north_m;  U = pe.up_m;
    H = hypot(E, N);
    t = pe.t_sec - min(pe.t_sec);

    hrmse = sqrt(mean(H.^2));   cep95 = prctile(H, 95);
    vrmse = sqrt(mean(U.^2));   vmean = mean(U);

    % ---- Figure 1: 수평/수직 오차 (2x2) ----------------------------------
    f1 = figure('Color', 'w', 'Name', 'WLS error vs truth', ...
                'Position', [80 80 1180 820]);

    subplot(2,2,1);
    scatter(E, N, 16, 'filled', 'MarkerFaceColor', RED, 'MarkerFaceAlpha', 0.5);
    hold on; grid on; axis equal; xline(0, 'k'); yline(0, 'k');
    xlabel('East error [m]'); ylabel('North error [m]');
    title(sprintf('Horizontal scatter  (RMSE %.2f m, CEP95 %.2f m)', hrmse, cep95));

    subplot(2,2,2);
    plot(t, H, '-', 'Color', RED, 'LineWidth', 0.8); hold on; grid on;
    yline(hrmse, '--k', sprintf('RMSE %.2f m', hrmse));
    xlabel('time since start [s]'); ylabel('horizontal error [m]');
    title('Horizontal error vs time');

    subplot(2,2,3);
    plot(t, U, '-', 'Color', BLUE, 'LineWidth', 0.8); hold on; grid on;
    yline(0, '-k'); yline(vmean, '--k', sprintf('mean %.1f m', vmean));
    xlabel('time since start [s]'); ylabel('vertical (up) error [m]');
    title(sprintf('Vertical error vs time  (RMSE %.1f m)', vrmse));

    subplot(2,2,4);
    Hs = sort(H);  Us = sort(abs(U));
    plot(Hs, linspace(0,100,numel(Hs)), '-', 'Color', RED, 'LineWidth', 1.4); hold on;
    plot(Us, linspace(0,100,numel(Us)), '-', 'Color', BLUE, 'LineWidth', 1.4);
    grid on; xlabel('error [m]'); ylabel('percentile [%]'); title('Error CDF');
    legend({'horizontal', 'vertical |U|'}, 'Location', 'southeast');

    sgtitle(sprintf('Positioning error vs truth  -  detector "%s"', best));
    exportgraphics(f1, fullfile(output_dir, 'matlab_error_horizontal_vertical.png'), ...
                   'Resolution', 130);

    % ---- Figure 2: 검출기 비교 (log) -------------------------------------
    [hr, order] = sort(hrmse_all);
    lbls = names(order);
    f2 = figure('Color', 'w', 'Name', 'Detector comparison', ...
                'Position', [80 80 960 500]);
    b = bar(hr, 'FaceColor', 'flat');
    for i = 1:numel(hr)
        b.CData(i,:) = (i == 1) .* RED + (i ~= 1) .* [0.5 0.55 0.6];
    end
    grid on; set(gca, 'YScale', 'log', 'XTick', 1:numel(hr), ...
        'XTickLabel', lbls, 'XTickLabelRotation', 22, 'TickLabelInterpreter', 'none');
    ylabel('Horizontal RMSE [m] (log)');
    title('WLS horizontal RMSE by outlier detector (vs truth)');
    ytop = max(hr) * 2.2;
    for i = 1:numel(hr)
        if hr(i) < 100, s = sprintf('%.2f', hr(i)); else, s = sprintf('%.0f km', hr(i)/1000); end
        text(i, hr(i), s, 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'FontSize', 8);
    end
    ylim([min(hr)*0.5, ytop]);
    exportgraphics(f2, fullfile(output_dir, 'matlab_detector_comparison.png'), ...
                   'Resolution', 130);

    % ---- Figure 3: 위성군 ablation (선택) --------------------------------
    abl_file = fullfile(output_dir, 'constellation_ablation.csv');
    if isfile(abl_file)
        abl = readtable(abl_file, 'VariableNamingRule', 'preserve');
        cfg = string(abl.config);  con = string(abl.constellation);
        tags = cfg + ":" + con;
        val = abl.horizontal_rmse_m;
        keep = ~isnan(val);
        f3 = figure('Color', 'w', 'Name', 'Constellation ablation', ...
                    'Position', [80 80 980 500]);
        bar(val(keep)); grid on; set(gca, 'YScale', 'log', ...
            'XTick', 1:nnz(keep), 'XTickLabel', tags(keep), 'XTickLabelRotation', 30, ...
            'TickLabelInterpreter', 'none');
        ylabel('Horizontal RMSE [m] (log)');
        title('Constellation ablation  (drop / only each constellation)');
        exportgraphics(f3, fullfile(output_dir, 'matlab_constellation_ablation.png'), ...
                       'Resolution', 130);
    end

    % ---- Figure 4: 측정 잡음 캘리브레이션 (선택) --------------------------
    cc_file = fullfile(output_dir, 'calibration_by_constellation.csv');
    if isfile(cc_file)
        cc = readtable(cc_file, 'VariableNamingRule', 'preserve');
        f4 = figure('Color', 'w', 'Name', 'Measurement noise calibration', ...
                    'Position', [80 80 1180 420]);

        subplot(1,3,1);
        bar(cc.clean_std_m); grid on; set(gca, 'YScale', 'log');
        set(gca, 'XTick', 1:height(cc), 'XTickLabel', string(cc.constellation), ...
            'XTickLabelRotation', 25, 'TickLabelInterpreter', 'none');
        ylabel('clean \sigma [m] (log)'); title('Per-constellation noise (vs truth)');
        for i = 1:height(cc)
            text(i, cc.clean_std_m(i), sprintf('%.0f', cc.clean_std_m(i)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
        end

        elf = fullfile(output_dir, 'calibration_by_elevation.csv');
        if isfile(elf)
            eb = readtable(elf, 'VariableNamingRule', 'preserve');
            xc = (eb.bin_lo + eb.bin_hi) / 2;
            subplot(1,3,2);
            plot(xc, eb.clean_std_m, '-o', 'Color', BLUE, 'LineWidth', 1.4, 'MarkerFaceColor', BLUE);
            grid on; xlabel('elevation [deg]'); ylabel('clean \sigma [m]');
            title('Noise vs elevation');
        end
        cnf = fullfile(output_dir, 'calibration_by_cn0.csv');
        if isfile(cnf)
            cb = readtable(cnf, 'VariableNamingRule', 'preserve');
            xc = (cb.bin_lo + cb.bin_hi) / 2;
            subplot(1,3,3);
            plot(xc, cb.clean_std_m, '-o', 'Color', RED, 'LineWidth', 1.4, 'MarkerFaceColor', RED);
            grid on; xlabel('C/N0 [dB-Hz]'); ylabel('clean \sigma [m]');
            title('Noise vs C/N0');
        end
        sgtitle('Empirical measurement noise (detrended residual vs truth)');
        exportgraphics(f4, fullfile(output_dir, 'matlab_calibration.png'), 'Resolution', 130);
    end

    % ---- Figure 5: sky plot + blunder catalog (선택) ---------------------
    cat_file = fullfile(output_dir, 'blunder_catalog.csv');
    if isfile(cat_file)
        cat = readtable(cat_file, 'VariableNamingRule', 'preserve');
        th = deg2rad(cat.azimuth_deg);
        rho = 90 - cat.elevation_deg;                 % 0=zenith(center), 90=horizon(edge)
        ar = abs(cat.residual_m);
        cap = prctile(ar, 95);
        isout = cat.is_outlier == 1;

        f5 = figure('Color', 'w', 'Name', 'Sky plot', 'Position', [80 80 820 720]);
        pax = polaraxes; hold(pax, 'on');
        pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
        polarscatter(pax, th, rho, 20, min(ar, cap), 'filled', 'MarkerFaceAlpha', 0.6);
        colormap(pax, turbo); cb = colorbar(pax);
        cb.Label.String = '|residual| [m] (capped at p95)';
        polarscatter(pax, th(isout), rho(isout), 70, 'r', 'x', 'LineWidth', 1.5);
        rlim(pax, [0 90]);
        pax.RTick = [0 30 60 90]; pax.RTickLabel = {'90', '60', '30', '0'};  % elevation
        title(pax, sprintf(['Sky plot  (N up, clockwise; ring = elevation)   ' ...
            '%d / %d flagged as outliers'], nnz(isout), height(cat)));
        exportgraphics(f5, fullfile(output_dir, 'matlab_skyplot.png'), 'Resolution', 130);

        fprintf('sky plot: %d/%d obs flagged as outliers\n', nnz(isout), height(cat));
    end

    % ---- Figure 6: ionosphere-treatment comparison (선택) ----------------
    ic_file = fullfile(output_dir, 'iono_comparison.csv');
    if isfile(ic_file)
        ic = readtable(ic_file, 'VariableNamingRule', 'preserve');
        modes = string(ic.mode);
        cols = lines(numel(modes));
        f6 = figure('Color', 'w', 'Name', 'Ionosphere comparison', 'Position', [60 60 1360 470]);

        subplot(1, 3, 1); hold on; grid on; axis equal;
        for i = 1:numel(modes)
            pe = fullfile(output_dir, "per_epoch_iono_" + modes(i) + ".csv");
            if ~isfile(pe), continue; end
            T = readtable(pe, 'VariableNamingRule', 'preserve');
            scatter(T.east_m, T.north_m, 8, cols(i, :), 'filled', 'MarkerFaceAlpha', 0.35);
        end
        xline(0, 'k'); yline(0, 'k');
        xlabel('East error [m]'); ylabel('North error [m]'); title('Horizontal error scatter');
        legend(modes, 'Location', 'bestoutside', 'Interpreter', 'none');

        subplot(1, 3, 2); hold on; grid on;
        for i = 1:numel(modes)
            pe = fullfile(output_dir, "per_epoch_iono_" + modes(i) + ".csv");
            if ~isfile(pe), continue; end
            T = readtable(pe, 'VariableNamingRule', 'preserve');
            hh = sort(hypot(T.east_m, T.north_m));
            plot(hh, linspace(0, 100, numel(hh)), 'Color', cols(i, :), 'LineWidth', 1.5);
        end
        xlabel('horizontal error [m]'); ylabel('percentile [%]'); title('Horizontal error CDF');
        legend(modes, 'Location', 'southeast', 'Interpreter', 'none');

        subplot(1, 3, 3);
        bar([ic.horizontal_rmse_m, ic.vertical_rmse_m]); grid on;
        set(gca, 'XTick', 1:numel(modes), 'XTickLabel', modes, 'XTickLabelRotation', 20, ...
            'TickLabelInterpreter', 'none');
        ylabel('RMSE [m]'); legend({'horizontal', 'vertical'}, 'Location', 'northwest');
        title('RMSE by ionosphere treatment');

        sgtitle('Horizontal nav solution vs ionosphere treatment');
        exportgraphics(f6, fullfile(output_dir, 'matlab_iono_comparison.png'), 'Resolution', 130);
    end

    fprintf('best detector: %s | H-RMSE %.2f m, CEP95 %.2f m, V-RMSE %.1f m\n', ...
            best, hrmse, cep95, vrmse);
    fprintf('saved matlab_*.png to %s\n', output_dir);
end
