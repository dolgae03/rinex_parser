base_station_path = "D:\과제\삼성 과제\Data\2차년도\Processed\BMHR20430089R_2025-09-28_20-03-48_processed.mat";
rover_path = "D:\과제\삼성 과제\Data\2차년도\Processed\29740_gnss_log_2025_09_28_17_00_44_processed.mat";

base_station = load(base_station_path);
rover = load(rover_path);

true_position = lla2ecef([36.372300, 127.358705, 90.594]);
% true_position = lla2ecef([36.372199, 127.3580734, 100.349579]);

true_velocity = [0, 0, 0];

%% geometric range 계산

tp = reshape(true_position, [1 1 3]);
% diff = base_station.XS_tot1 - tp;     % [Nt x Ns x 3]
diff = base_station.XS_tot1_rot - tp;
ranges = sqrt(sum(diff.^2, 3));       % [Nt x Ns]
% 보정 정보c
base_station.pr_correction = ranges - base_station.pr1;

%% doppler correction 계산
los = diff ./ ranges; 

rv = reshape(true_velocity, [1 1 3]);                         % [1 x 1 x 3]
rel_vel = base_station.VS_tot1_rot - rv;

C = 299792458;  % [m/s]
lambda = C ./ 1575.42e6;   

range_rate = sum(rel_vel .* los, 3);

% 보정 정보
geom_dopp_hz = - range_rate ./ lambda;      
base_station.dop_correction = geom_dopp_hz - base_station.dop1;  % [Nt x Ns]  
% [Nt x Ns], [Hz]
mean_after = mean(geom_dopp_hz(:) - (base_station.dop1(:) + base_station.dop_correction(:)), 'omitnan');

%% smartphone data to station
[rover, base_station] = timesync(rover, base_station);
corrected_rover = measurement_copy(rover, base_station);

%% 저장 경로 만들기 (_correction_added 붙이기)
[pathstr, name, ext] = fileparts(rover_path);
save_path = fullfile(pathstr, name + "_correction_added" + ext);

% 저장
save(save_path, "-struct", "corrected_rover");