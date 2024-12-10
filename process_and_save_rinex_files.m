function process_and_save_rinex_files(input_dir, nav_dir, output_dir)
    % Suppress all warnings
    warning('off', 'all');
    
    % Add paths for necessary directories
    addpath(genpath(input_dir));

    file_extension = '.24o';  % 찾을 파일의 확장자
    rinex_files = find_all_files_with_extension(input_dir, file_extension);
    
    % Process each file and save output
    for i = 1:length(rinex_files)
        file_path = rinex_files{i};
        
        % Extract the relative path from input_dir to the file
        relative_path = strrep(file_path, [input_dir filesep], '');
        [relative_folder, file_name, ~] = fileparts(relative_path);
        
        % Create the corresponding output folder path
        output_folder_path = fullfile(output_dir, relative_folder);
        
        % Make sure the output folder exists
        if ~exist(output_folder_path, 'dir')
            mkdir(output_folder_path);
        end
        
        % Define the output file path with the original folder structure
        output_file_path = fullfile(output_folder_path, [file_name, '_processed.mat']);
        
        % Process the RINEX observation data
        data = process_rinex_file(file_path, nav_dir);
        
        % Save the processed data
        save(output_file_path, '-struct', 'data');
    end
end

function rinex_files = find_all_files_with_extension(root_dir, extension)
    % 특정 디렉토리와 모든 하위 디렉토리에서 주어진 확장자의 파일을 찾음
    % root_dir: 탐색을 시작할 최상위 디렉토리
    % extension: 찾을 파일의 확장자 (예: '.21o')
    
    % 초기화
    rinex_files = {};  % 셀 배열로 초기화하여 각 파일 경로를 저장
    
    % 현재 디렉토리의 모든 파일과 폴더 가져오기
    items = dir(root_dir);
    
    % 각 항목 검사
    for i = 1:length(items)
        % 현재 항목의 전체 경로
        current_path = fullfile(root_dir, items(i).name);
        
        % 디렉토리인지 검사
        if items(i).isdir
            % '.' 와 '..' 디렉토리는 건너뜀
            if ~ismember(items(i).name, {'.', '..'})
                % 하위 디렉토리에서 재귀적으로 호출
                subdir_files = find_all_files_with_extension(current_path, extension);
                rinex_files = [rinex_files, subdir_files];  % 셀 배열 결합
            end
        elseif endsWith(items(i).name, extension)
            % 파일 확장자가 일치하는 경우, 파일 경로 추가
            rinex_files{end + 1} = current_path;
        end
    end
end


function data = process_rinex_file(file_path, nav_file_path)
    % Set up directories and add necessary paths
    homedir = pwd;
    source_dir = sprintf('%s/source', homedir);
    addpath(genpath(source_dir));
    custom_dir = sprintf('%s/custom', homedir);
    addpath(genpath(custom_dir));
    navutils_dir = sprintf('%s/navutils', homedir);
    addpath(genpath(navutils_dir));

    % Load the RINEX observation data
    cc = Constellation_Collector([1, 1, 1, 1, 1, 0, 0]);

    processing_interval = 1;  % Adjust as needed
    [pr1, ph1, pr2, ph2, pr3, ph3, dop1, dop2, dop3, snr1, snr2, snr3, ...
        time_zero, time_GPS, time, week, date, pos, interval, antoff, ...
        antmod, codeC1, marker] = load_RINEX_obs(file_path, cc, processing_interval);

    % Initialize arrays for satellite positions and velocities
    nSatTot = cc.getNumSat();
    XS_tot1 = NaN(size(pr1, 2), nSatTot, 3);
    VS_tot1 = NaN(size(pr1, 2), nSatTot, 3);

    %LEAP SECONDS for gps time
    % LEAP_SECOND = 18
    % time_zero = time_zero + LEAP_SECOND;

    % Set up other parameters
    nsat = size(pr1, 1);
    time_rx = time_GPS + time_zero;
    dtR = 0;
    err_tropo = zeros(nsat, 1);
    err_iono = zeros(nsat, 1);
    p_rate = 1e-6;
    

    % Load the current navigation file
    [Eph, iono] = load_RINEX_nav(nav_file_path, cc, 0, 2);
    lambda = goGNSS.getGNSSWavelengths(Eph, [], nSatTot);

    % Calculate satellite positions and velocities for each time step
    for i = 1:size(pr1, 2)
        avail_sat1 = find(pr1(:, i)); % Find available satellites for the time step
        [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, eclipsed, sys_idx] = ...
            satellite_positions(time_rx(i), pr1(:, i), avail_sat1, Eph, [], [], ...
            err_tropo, err_iono, dtR, 1, 'NONE', lambda, p_rate);

        for j = 1:size(XS_tx, 1)
            % Update XS_tot1 and VS_tot1 with positions and velocities
            XS_tot1(i, avail_sat1(j), :) = XS_tx(j, :); % at transmission time
            VS_tot1(i, avail_sat1(j), :) = VS_tx(j, :);
        end
    end

    [pr1, ph1, pr2, ph2, pr3, ph3, dop1, dop2, dop3, snr1, snr2, snr3, XS_tot1, VS_tot1] = ...
          zero_nan(pr1, ph1, pr2, ph2, pr3, ph3, dop1, dop2, dop3, snr1, snr2, snr3, XS_tot1, VS_tot1);

    % Organize data into a struct for saving
    data.pr1 = pr1.';
    data.ph1 = ph1.';
    data.pr2 = pr2.';
    data.ph2 = ph2.';
    data.pr3 = pr3.';
    data.ph3 = ph3.';
    data.dop1 = dop1.';
    data.dop2 = dop2.';
    data.dop3 = dop3.';
    data.snr1 = snr1.';
    data.snr2 = snr2.';
    data.snr3 = snr3.';
    data.time_zero = time_zero;
    data.time_GPS = time_GPS;
    data.time = time;
    data.week = week;
    data.date = date;
    data.pos = pos;
    data.interval = interval;
    data.antoff = antoff;
    data.antmod = antmod;
    data.codeC1 = codeC1;
    data.marker = marker;
    data.active_constellation = cc.active_list;
    data.constellation_name = cc.SYS_NAME;
    data.constellation_idx = cc.IDX_SAT;

    % Add satellite positions and velocities to the data structure
    data.XS_tot1 = XS_tot1;
    data.VS_tot1 = VS_tot1;

    data.SVpos_x = squeeze(XS_tot1(:,:,1)); data.SVpos_y = squeeze(XS_tot1(:,:,2)); data.SVpos_z = squeeze(XS_tot1(:,:,3));
    data.SVvel_x = squeeze(VS_tot1(:,:,1)); data.SVvel_y = squeeze(VS_tot1(:,:,2)); data.SVvel_z = squeeze(VS_tot1(:,:,3));
end
