% Define the file path
% filePath = "D:\과제\삼성 과제\Data\Processed\1차년도 Data 수집\reciever_opensky_2025_04_09_processed.mat";
% filePath = "D:\과제\삼성 과제\Data\2차년도\Processed\BMHR21030714D_2025-09-19_15-25-44_processed.mat";
% filePath = "D:\과제\삼성 과제\Data\2차년도\Processed\receiver_opensky_0514_processed.mat";
% filePath = "D:\과제\삼성 과제\Data\2차년도\Processed\BMHR20430089R_2025-09-28_20-03-48_processed_correction_added.mat";
filePath = "D:\과제\삼성 과제\Data\2차년도\Processed\29740_gnss_log_2025_09_28_17_00_44_processed_correction_added.mat";
% filePath = "D:\과제\삼성 과제\Data\2차년도\Processed\Garden_Data\BMHR20430089R_2025-10-01_14-23-09_processed_correction_added.mat";
% filePath = "D:\과제\삼성 과제\Data\2차년도\Processed\SPP_Check_2025-10-10_12-08-11_processed.mat";
% filePath = "D:\과제\삼성 과제\Data\2차년도\Processed\BMHR21030714D_2025-10-11_10-57-32_processed.mat";

% Load the .mat file
data = load(filePath);

[pathstr, name, ext] = fileparts(filePath);
save_path = name + '.csv';

export_gnss_to_csv_fast(filePath, "./data/" + save_path)

function export_gnss_to_csv_fast(mat_file, out_csv)
    data = load(mat_file);

    % 시간 관련
    time_base = data.time_zero;
    time_sec  = data.time_GPS(:);
    gps_week  = data.week(:);
    tow_sec   = time_base + time_sec - gps_week*7*24*3600;

    n_epoch = numel(time_sec);
    n_sv    = size(data.SVpos_x_rot, 2);
    n_freq  = 3;

    % --- 주파수별 관측데이터 쌓기 ---
    pr_all  = cat(3, data.pr1, data.pr2, data.pr3);
    ph_all  = cat(3, data.ph1, data.ph2, data.ph3);
    dop_all = cat(3, data.dop1, data.dop2, data.dop3);
    snr_all = cat(3, data.snr1, data.snr2, data.snr3);

    if isfield(data, "pr_correction")
        pr_corr = data.pr_correction;
    else
        pr_corr = zeros(size(data.pr1));
    end
    
    if isfield(data, "dop_correction")
        dop_corr = data.dop_correction;
    else
        dop_corr = zeros(size(data.pr1));
    end

    loi = zeros(size(data.pr1));

    % --- 전체 인덱스 생성 ---
    [I, J, K] = ndgrid(1:n_epoch, 1:n_sv, 1:n_freq);
    I = I(:); J = J(:); K = K(:);

    pr_val = pr_all(sub2ind(size(pr_all), I, J, K));

    valid_mask = ~isnan(pr_val);
    I = I(valid_mask); J = J(valid_mask); K = K(valid_mask);
    pr_val = pr_val(valid_mask);

    N = numel(I);

    % --- 별자리 매핑 ---
    edges = [data.constellation_idx(:).', n_sv+1];
    mapping_to_format = [0 4 1 3 2 5 6];

    constellation = zeros(N,1);
    prn = zeros(N,1);
    for g = 1:numel(mapping_to_format)
        mask = J >= edges(g) & J < edges(g+1);
        constellation(mask) = mapping_to_format(g);
        prn(mask) = J(mask) - edges(g) + 1;
    end

    % --- 테이블 컬럼 벡터 생성 ---
    T = table;
    T.t_sec    = time_base + time_sec(I);
    T.gps_week = gps_week(I);
    T.tow_sec  = tow_sec(I);
    T.constellation = constellation;
    T.prn = prn;
    T.sv_pos_x = data.SVpos_x(sub2ind([n_epoch,n_sv], I, J));
    T.sv_pos_y = data.SVpos_y(sub2ind([n_epoch,n_sv], I, J));
    T.sv_pos_z = data.SVpos_z(sub2ind([n_epoch,n_sv], I, J));
    T.sv_vel_x = data.SVvel_x(sub2ind([n_epoch,n_sv], I, J));
    T.sv_vel_y = data.SVvel_y(sub2ind([n_epoch,n_sv], I, J));
    T.sv_vel_z = data.SVvel_z(sub2ind([n_epoch,n_sv], I, J));
    T.sv_clock_bias = data.sv_clock_bias(sub2ind([n_epoch,n_sv], I, J));
    T.sv_clock_drift = zeros(N,1);
    T.pr_correction = pr_corr(sub2ind([n_epoch,n_sv], I, J));
    T.dop_correction = dop_corr(sub2ind([n_epoch,n_sv], I, J));
    freq_list = [1575.42e6, 1227.60e6, 1176.45e6];
    T.frequency_hz = freq_list(K).';   % ← 전치(.')로 N×1 맞추기
    T.code_type = repmat("C", N, 1);
    T.pseudorange_m = pr_val;
    T.phase_cycle   = ph_all(sub2ind(size(pr_all), I, J, K));
    T.doppler_hz    = dop_all(sub2ind(size(dop_all), I, J, K));
    T.snr_dbhz      = snr_all(sub2ind(size(snr_all), I, J, K));
    T.loi = loi(sub2ind([n_epoch,n_sv], I, J));

    iono = num2cell(repmat(data.iono(:).', N, 1), 2);
    iono_mat = cell2mat(iono);
    T.iono_a0 = iono_mat(:,1); T.iono_a1 = iono_mat(:,2);
    T.iono_a2 = iono_mat(:,3); T.iono_a3 = iono_mat(:,4);
    T.iono_b0 = iono_mat(:,5); T.iono_b1 = iono_mat(:,6);
    T.iono_b2 = iono_mat(:,7); T.iono_b3 = iono_mat(:,8);

    % --- 🔹 시간순 정렬 (t_sec 기준 오름차순) ---
    T = sortrows(T, ["t_sec", "constellation", "prn"]);

    varNames = {'t_sec','gps_week','tow_sec','constellation','prn', ...
                 'sv_pos_x','sv_pos_y','sv_pos_z', ...
                 'sv_vel_x','sv_vel_y','sv_vel_z', ...
                 'iono_a0','iono_a1','iono_a2','iono_a3', ...
                 'iono_b0','iono_b1','iono_b2','iono_b3', ...
                 'sv_clock_bias','sv_clock_drift', ...
                 'pr_correction', 'dop_correction', ...
                 'frequency_hz','code_type', ...
                 'pseudorange_m','phase_cycle','doppler_hz','snr_dbhz','loi'};
    
    T = T(:, varNames);  % ✅ 순서 강제
    % --- 저장 ---
    writetable(T, out_csv, "Delimiter", '\t');
    fprintf("CSV saved to %s (rows=%d)\n", out_csv, height(T));
end

