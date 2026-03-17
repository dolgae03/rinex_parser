function summary = fill_tsv_with_rinex_nav(tsv_file, nav_file, out_tsv)
    validateattributes(tsv_file, {'char', 'string'}, {'scalartext'});
    validateattributes(nav_file, {'char', 'string'}, {'scalartext'});
    validateattributes(out_tsv, {'char', 'string'}, {'scalartext'});

    tsv_file = char(tsv_file);
    nav_file = char(nav_file);
    out_tsv = char(out_tsv);

    if ~exist(tsv_file, 'file')
        error('TSV file not found: %s', tsv_file);
    end
    if ~exist(nav_file, 'file')
        error('Navigation RINEX file not found: %s', nav_file);
    end

    opts = detectImportOptions(tsv_file, 'FileType', 'text', 'Delimiter', '\t');
    opts.VariableNamingRule = 'preserve';
    T = readtable(tsv_file, opts);

    required_vars = {'t_sec', 'gps_week', 'tow_sec', 'constellation', 'prn', 'pseudorange_m'};
    missing_vars = required_vars(~ismember(required_vars, T.Properties.VariableNames));
    if ~isempty(missing_vars)
        error('Missing required TSV columns: %s', strjoin(missing_vars, ', '));
    end

    target_vars = {'sv_pos_x', 'sv_pos_y', 'sv_pos_z', ...
                   'sv_vel_x', 'sv_vel_y', 'sv_vel_z', ...
                   'sv_clock_bias', 'sv_clock_drift', ...
                   'iono_a0', 'iono_a1', 'iono_a2', 'iono_a3', ...
                   'iono_b0', 'iono_b1', 'iono_b2', 'iono_b3'};
    for i = 1:numel(target_vars)
        if ~ismember(target_vars{i}, T.Properties.VariableNames)
            T.(target_vars{i}) = nan(height(T), 1);
        end
    end

    % Keep goGPS's historical fixed satellite indexing layout.
    cc = Constellation_Collector([1, 1, 1, 1, 1, 0, 0]);
    [Eph, iono] = load_RINEX_nav(nav_file, cc, 0, 0);
    if isempty(Eph) || size(Eph, 1) < 33
        error(['Navigation parsing did not produce a valid ephemeris matrix for %s. ' ...
               'For mixed navigation files, goGPS expects a filename ending in ''p'' (for example .25p).'], nav_file);
    end
    T = fill_iono_coefficients(T, iono);

    nSatTot = cc.getNumSat();
    lambda = goGNSS.getGNSSWavelengths(Eph, [], nSatTot);

    epoch_keys = double(T.t_sec);
    [unique_epochs, ~, epoch_groups] = unique(epoch_keys, 'stable');

    for epoch_idx = 1:numel(unique_epochs)
        row_idx = find(epoch_groups == epoch_idx);
        time_rx = unique_epochs(epoch_idx);

        [sat_list, sat_row_map, pseudorange_vec] = build_epoch_satellite_inputs(T, row_idx, cc, nSatTot);
        if isempty(sat_list)
            continue;
        end

        err_tropo = zeros(nSatTot, 1);
        err_iono = zeros(nSatTot, 1);
        dtR = 0;
        p_rate = 1e-6;

        [~, ~, dtS, XS_tx, VS_tx, time_tx, no_eph, ~, ~, eph_rows] = ...
            satellite_positions(time_rx, pseudorange_vec, sat_list, Eph, [], [], ...
                                err_tropo, err_iono, dtR, 1, 'NONE', lambda, p_rate);

        for j = 1:numel(sat_list)
            sat_id = sat_list(j);
            if no_eph(j)
                continue;
            end

            rows_for_sat = sat_row_map{sat_id};
            if isempty(rows_for_sat)
                continue;
            end

            T = fill_missing_scalar(T, 'sv_pos_x', rows_for_sat, XS_tx(j, 1));
            T = fill_missing_scalar(T, 'sv_pos_y', rows_for_sat, XS_tx(j, 2));
            T = fill_missing_scalar(T, 'sv_pos_z', rows_for_sat, XS_tx(j, 3));
            T = fill_missing_scalar(T, 'sv_vel_x', rows_for_sat, VS_tx(j, 1));
            T = fill_missing_scalar(T, 'sv_vel_y', rows_for_sat, VS_tx(j, 2));
            T = fill_missing_scalar(T, 'sv_vel_z', rows_for_sat, VS_tx(j, 3));
            T = fill_missing_scalar(T, 'sv_clock_bias', rows_for_sat, dtS(j) * goGNSS.V_LIGHT);
            T = fill_missing_scalar(T, 'sv_clock_drift', rows_for_sat, ...
                compute_satellite_clock_drift(time_tx(j), eph_rows{j}) * goGNSS.V_LIGHT);
        end
    end

    out_dir = fileparts(out_tsv);
    if ~isempty(out_dir) && ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    writetable(T, out_tsv, 'Delimiter', '\t', 'FileType', 'text');

    summary.output_tsv = out_tsv;
    summary.nav_file = nav_file;
    summary.total_rows = height(T);
    summary.finite_sv_pos_rows = nnz(isfinite(T.sv_pos_x) & isfinite(T.sv_pos_y) & isfinite(T.sv_pos_z));
    summary.finite_sv_vel_rows = nnz(isfinite(T.sv_vel_x) & isfinite(T.sv_vel_y) & isfinite(T.sv_vel_z));
    summary.finite_sv_clock_bias_rows = nnz(isfinite(T.sv_clock_bias));
    summary.finite_sv_clock_drift_rows = nnz(isfinite(T.sv_clock_drift));
    summary.filled_iono_rows = nnz(isfinite(T.iono_a0) | isfinite(T.iono_b0));
end

function [sat_list, sat_row_map, pseudorange_vec] = build_epoch_satellite_inputs(T, row_idx, cc, nSatTot)
    pseudorange_vec = nan(nSatTot, 1);
    sat_row_map = cell(nSatTot, 1);

    const_values = double(T.constellation(row_idx));
    prn_values = double(T.prn(row_idx));
    pr_values = double(T.pseudorange_m(row_idx));

    if ismember('snr_dbhz', T.Properties.VariableNames)
        snr_values = double(T.snr_dbhz(row_idx));
    else
        snr_values = nan(numel(row_idx), 1);
    end

    keys = [const_values, prn_values];
    [unique_keys, ~, key_groups] = unique(keys, 'rows', 'stable');

    for key_idx = 1:size(unique_keys, 1)
        sub_rows = row_idx(key_groups == key_idx);
        constellation = unique_keys(key_idx, 1);
        prn = unique_keys(key_idx, 2);
        sat_id = map_to_gogps_satellite_index(constellation, prn, cc);
        if isnan(sat_id) || sat_id < 1 || sat_id > nSatTot
            continue;
        end

        local_mask = key_groups == key_idx;
        local_pr = pr_values(local_mask);
        local_snr = snr_values(local_mask);
        candidate_local_idx = find(isfinite(local_pr) & local_pr > 0);
        if isempty(candidate_local_idx)
            continue;
        end

        if any(isfinite(local_snr(candidate_local_idx)))
            [~, best_rel] = max(local_snr(candidate_local_idx));
            best_local = candidate_local_idx(best_rel);
        else
            best_local = candidate_local_idx(1);
        end

        representative_pr = local_pr(best_local);
        pseudorange_vec(sat_id) = representative_pr;
        sat_row_map{sat_id} = sub_rows;
    end

    sat_list = find(isfinite(pseudorange_vec));
end

function sat_id = map_to_gogps_satellite_index(constellation, prn, cc)
    switch constellation
        case 0
            sat_id = cc.IDX_SAT(cc.ID_GPS) + prn - 1;
        case 1
            sat_id = cc.IDX_SAT(cc.ID_GALILEO) + prn - 1;
        case 2
            sat_id = cc.IDX_SAT(cc.ID_BEIDOU) + prn - 1;
        case 3
            sat_id = cc.IDX_SAT(cc.ID_GLONASS) + prn - 1;
        case 4
            sat_id = cc.IDX_SAT(cc.ID_QZSS) + prn - 1;
        otherwise
            sat_id = nan;
    end
end

function T = fill_iono_coefficients(T, iono)
    if isempty(iono) || numel(iono) < 8 || ~any(isfinite(double(iono)))
        return;
    end

    coeff_names = {'iono_a0', 'iono_a1', 'iono_a2', 'iono_a3', ...
                   'iono_b0', 'iono_b1', 'iono_b2', 'iono_b3'};
    coeff_values = double(iono(:)');
    coeff_values = coeff_values(1:8);

    for idx = 1:numel(coeff_names)
        current_values = double(T.(coeff_names{idx}));
        replace_mask = ~isfinite(current_values) | current_values == 0;
        if any(replace_mask) && isfinite(coeff_values(idx))
            T.(coeff_names{idx})(replace_mask) = coeff_values(idx);
        end
    end
end

function T = fill_missing_scalar(T, var_name, rows, value)
    if ~isfinite(value)
        return;
    end

    current_values = double(T.(var_name)(rows));
    replace_mask = ~isfinite(current_values);
    if any(replace_mask)
        target_rows = rows(replace_mask);
        T.(var_name)(target_rows) = value;
    end
end

function drift = compute_satellite_clock_drift(time_value, eph_row)
    if isempty(eph_row) || any(~isfinite(eph_row))
        drift = nan;
        return;
    end

    sys_id = char(eph_row(31));
    if sys_id == 'R'
        drift = eph_row(3);
        return;
    end

    af2 = eph_row(2);
    af1 = eph_row(20);
    ref_toc = eph_row(33);

    if sys_id == 'C'
        time_value = time_value - 14;
    end

    dt = check_t(time_value - ref_toc);
    drift = af1 + 2 * af2 * dt;
end
