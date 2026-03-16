function summary = export_processed_mat_to_tsv(mat_file, out_csv, display_rows)
    if nargin < 3 || isempty(display_rows)
        display_rows = 10;
    end

    data = load(mat_file);

    required_fields = {'time_zero', 'time_GPS', 'week', 'pr1', 'SVpos_x', 'SVpos_y', 'SVpos_z', ...
                       'SVvel_x', 'SVvel_y', 'SVvel_z', 'sv_clock_bias', 'constellation_idx'};
    for i = 1:numel(required_fields)
        if ~isfield(data, required_fields{i})
            error('Missing required field "%s" in %s', required_fields{i}, mat_file);
        end
    end

    time_base = data.time_zero;
    time_sec = data.time_GPS(:);
    gps_week = data.week(:);
    tow_sec = time_base + time_sec - gps_week * 7 * 24 * 3600;

    n_epoch = numel(time_sec);
    n_sv = size(data.SVpos_x, 2);
    n_freq = 3;

    pr1 = double(data.pr1);
    pr2 = get_or_nan(data, 'pr2', size(pr1));
    pr3 = get_or_nan(data, 'pr3', size(pr1));
    ph1 = get_or_nan(data, 'ph1', size(pr1));
    ph2 = get_or_nan(data, 'ph2', size(pr1));
    ph3 = get_or_nan(data, 'ph3', size(pr1));
    dop1 = get_or_nan(data, 'dop1', size(pr1));
    dop2 = get_or_nan(data, 'dop2', size(pr1));
    dop3 = get_or_nan(data, 'dop3', size(pr1));
    snr1 = get_or_nan(data, 'snr1', size(pr1));
    snr2 = get_or_nan(data, 'snr2', size(pr1));
    snr3 = get_or_nan(data, 'snr3', size(pr1));

    pr_all = cat(3, pr1, pr2, pr3);
    ph_all = cat(3, ph1, ph2, ph3);
    dop_all = cat(3, dop1, dop2, dop3);
    snr_all = cat(3, snr1, snr2, snr3);

    pr_corr = get_or_nan(data, 'pr_correction', size(pr1));
    has_dop_corr = isfield(data, 'dop_correction');
    if has_dop_corr
        dop_corr = get_or_nan(data, 'dop_correction', size(pr1));
    end

    loi = zeros(size(pr1));

    [I, J, K] = ndgrid(1:n_epoch, 1:n_sv, 1:n_freq);
    I = I(:);
    J = J(:);
    K = K(:);

    pr_val = pr_all(sub2ind(size(pr_all), I, J, K));
    valid_mask = ~isnan(pr_val);

    I = I(valid_mask);
    J = J(valid_mask);
    K = K(valid_mask);
    pr_val = pr_val(valid_mask);

    N = numel(I);

    edges = [data.constellation_idx(:).', n_sv + 1];
    mapping_to_format = [0 4 1 3 2 5 6];

    constellation = zeros(N, 1);
    prn = zeros(N, 1);
    for g = 1:numel(mapping_to_format)
        mask = J >= edges(g) & J < edges(g + 1);
        constellation(mask) = mapping_to_format(g);
        prn(mask) = J(mask) - edges(g) + 1;
    end

    T = table;
    T.t_sec = time_base + time_sec(I);
    T.gps_week = gps_week(I);
    T.tow_sec = tow_sec(I);
    T.constellation = constellation;
    T.prn = prn;
    T.sv_pos_x = data.SVpos_x(sub2ind([n_epoch, n_sv], I, J));
    T.sv_pos_y = data.SVpos_y(sub2ind([n_epoch, n_sv], I, J));
    T.sv_pos_z = data.SVpos_z(sub2ind([n_epoch, n_sv], I, J));
    T.sv_vel_x = data.SVvel_x(sub2ind([n_epoch, n_sv], I, J));
    T.sv_vel_y = data.SVvel_y(sub2ind([n_epoch, n_sv], I, J));
    T.sv_vel_z = data.SVvel_z(sub2ind([n_epoch, n_sv], I, J));
    [iono_a0, iono_a1, iono_a2, iono_a3, iono_b0, iono_b1, iono_b2, iono_b3] = build_iono_columns(data, N);
    T.iono_a0 = iono_a0;
    T.iono_a1 = iono_a1;
    T.iono_a2 = iono_a2;
    T.iono_a3 = iono_a3;
    T.iono_b0 = iono_b0;
    T.iono_b1 = iono_b1;
    T.iono_b2 = iono_b2;
    T.iono_b3 = iono_b3;
    T.sv_clock_bias = data.sv_clock_bias(sub2ind([n_epoch, n_sv], I, J));
    T.sv_clock_drift = nan(N, 1);
    T.pr_correction = pr_corr(sub2ind([n_epoch, n_sv], I, J));
    if has_dop_corr
        T.dop_correction = dop_corr(sub2ind([n_epoch, n_sv], I, J));
    end
    freq_list = [1575.42e6, 1227.60e6, 1176.45e6];
    T.frequency_hz = freq_list(K).';
    T.code_type = repmat("C", N, 1);
    T.pseudorange_m = pr_val;
    T.phase_cycle = ph_all(sub2ind(size(ph_all), I, J, K));
    T.doppler_hz = dop_all(sub2ind(size(dop_all), I, J, K));
    T.snr_dbhz = snr_all(sub2ind(size(snr_all), I, J, K));
    T.loi = loi(sub2ind([n_epoch, n_sv], I, J));
    T.gt_pos_x = nan(N, 1);
    T.gt_pos_y = nan(N, 1);
    T.gt_pos_z = nan(N, 1);

    T = sortrows(T, {'t_sec', 'constellation', 'prn', 'frequency_hz'});
    writetable(T, out_csv, 'Delimiter', '\t');

    finite_sv_pos_mask = isfinite(T.sv_pos_x) & isfinite(T.sv_pos_y) & isfinite(T.sv_pos_z);
    finite_sv_vel_mask = isfinite(T.sv_vel_x) & isfinite(T.sv_vel_y) & isfinite(T.sv_vel_z);
    finite_sv_clock_mask = isfinite(T.sv_clock_bias);

    preview_rows = min(display_rows, height(T));
    if preview_rows > 0
        summary.preview = T(1:preview_rows, {'t_sec', 'constellation', 'prn', ...
                                             'sv_pos_x', 'sv_pos_y', 'sv_pos_z', ...
                                             'sv_vel_x', 'sv_vel_y', 'sv_vel_z'});
    else
        summary.preview = T([], :);
    end

    summary.mat_file = mat_file;
    summary.out_csv = out_csv;
    summary.row_count = height(T);
    summary.epoch_count = n_epoch;
    summary.finite_sv_pos_rows = nnz(finite_sv_pos_mask);
    summary.finite_sv_vel_rows = nnz(finite_sv_vel_mask);
    summary.finite_sv_clock_rows = nnz(finite_sv_clock_mask);
end

function value = get_or_nan(data, field_name, fallback_size)
    if isfield(data, field_name)
        value = double(data.(field_name));
    else
        value = nan(fallback_size);
    end
end

function [a0, a1, a2, a3, b0, b1, b2, b3] = build_iono_columns(data, row_count)
    if isfield(data, 'iono') && numel(data.iono) >= 8
        iono = double(data.iono(:));
        iono = iono(1:8).';
        iono_mat = repmat(iono, row_count, 1);
    else
        iono_mat = nan(row_count, 8);
    end

    a0 = iono_mat(:, 1);
    a1 = iono_mat(:, 2);
    a2 = iono_mat(:, 3);
    a3 = iono_mat(:, 4);
    b0 = iono_mat(:, 5);
    b1 = iono_mat(:, 6);
    b2 = iono_mat(:, 7);
    b3 = iono_mat(:, 8);
end
