%% 입력 (이미 가지고 있는 값들)
parsed_lla_pos = parse_bestposa_lla("C:\Users\mskim\Documents\NovAtelData\DataLogger\BMHR21030714D_2025-10-11_10-57-32\Converted_On_20251011_T1426\ASCII\BMHR21030714D_2025-10-11_10-57-32_BESTPOS.ASCII");

addpath("navutils\");   % wgsxyz2lla, wgslla2enu 등

true_pose_ecef = [-3119951.884, 4086922.291, 3761584.638];   % 기준점 (ECEF, m)
[ref_lat, ref_lon, ref_alt] = wgsxyz2lla(true_pose_ecef);    % 기준점 LLA

%% 설정: 지오이드 보정 (한국 기준 대략 +30 m)
GEOID_OFFSET_M = 25.0;   % WGS84 타원체고에 +30 m 더해서 MSL 근사
% 만약 BESTPOS altitude가 이미 HAE(타원체)이고, ENU 기준도 타원체 기준으로 하고 싶다면 0으로 바꾸세요.
% GEOID_OFFSET_M = 0;

%% ENU 변환 (벡터화)
n = size(parsed_lla_pos, 1);
enu = nan(n, 3);   % [E N U]

% parsed_lla_pos(:,1:3) = [lat deg, lon deg, alt m(HAE)]
lat_vec = parsed_lla_pos(:,1);
lon_vec = parsed_lla_pos(:,2);
alt_vec = parsed_lla_pos(:,3) + GEOID_OFFSET_M;  % 지오이드 보정 적용

% wgslla2enu가 스칼라용이면 loop/arrayfun 사용
for i = 1:n
    enu(i,:) = wgslla2enu(lat_vec(i), lon_vec(i), alt_vec(i), ...
                          ref_lat, ref_lon, ref_alt);
end

E = enu(:,1); N = enu(:,2); U = enu(:,3);
R2 = sqrt(E.^2 + N.^2);      % 2D 수평 거리
R3 = sqrt(E.^2 + N.^2 + U.^2);

%% 요약 통계
meanE = mean(E, 'omitnan');     
stdE  = std(E, 1, 'omitnan');     
rmsE  = sqrt(mean(E.^2, 'omitnan'));

meanN = mean(N, 'omitnan');     
stdN  = std(N, 1, 'omitnan');     
rmsN  = sqrt(mean(N.^2, 'omitnan'));

meanU = mean(U, 'omitnan');     
stdU  = std(U, 1, 'omitnan');     
rmsU  = sqrt(mean(U.^2, 'omitnan'));

rms2D = mean(R2, 'omitnan');
d2rms = 2*sqrt(var(E,1,'omitnan')+var(N,1,'omitnan'));

rms3D = mean(R3, 'omitnan');
d3rms = mean(R3);

% CEP(경험적 분위수)
cep50 = prctile(R2, 50);
cep95 = prctile(R2, 95);

%% 그림 1: EN 평면 산점도
figure('Name', 'EN Scatter (Horizontal Error)', 'Color', 'w');
scatter(E, N, 12, 'filled', 'MarkerFaceColor','r'); hold on; grid on; axis equal;
plot(meanE, meanN, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y'); % 평균점
xlabel('East [m]'); ylabel('North [m]');
title('Horizontal Scatter (E vs N)');

% % 2DRMS 원 그리기
% th = linspace(0, 2*pi, 360);
% plot(d2rms/2*cos(th), d2rms/2*sin(th), 'r-', 'LineWidth', 1.2); % 반지름 = d2rms/2? → 주의
% % 참고: 2DRMS 정의가 표준에 따라 다를 수 있음. 일반적으로 2*sqrt(sigmaE^2+sigmaN^2)가 “2DRMS”이며
% % 그 값 자체가 지름이 아니라 반지름으로 쓰는 경우도 많아 혼동됨.
% % 여기서는 "2DRMS 값을 반지름으로 갖는 원"을 그리려면 아래가 맞음:
% plot(d2rms*cos(th), d2rms*sin(th), 'r-', 'LineWidth', 1.2); % 반지름 = 2DRMS

legend('Samples', 'Mean', '2DRMS circle', 'Location', 'bestoutside');

txt = sprintf(['Mean [E N]=[%.3f %.3f] m\n',...
               'STD  [E N]=[%.3f %.3f] m\n',...
               'RMS  [E N]=[%.3f %.3f] m\n',...
               'RMS2D=%.3f m,  2DRMS=%.3f m\n',...
               'CEP50=%.3f m, CEP95=%.3f m'], ...
               meanE, meanN, stdE, stdN, rmsE, rmsN, rms2D, d2rms, cep50, cep95);
annotation('textbox',[0.72 0.60 0.25 0.3],'String',txt,'FitBoxToText','on','BackgroundColor','w');

%% 그림 2: 시간축 E/N/U 오차
t = (1:n).';   % 시간 인덱스 (초단위 timestamp 없으니 시퀀스 인덱스 사용)

figure('Name', 'ENU Error vs Epoch', 'Color', 'w');

subplot(3,1,1);
plot(t, E, 'LineWidth', 1); grid on; ylabel('E [m]'); title('E Error');
yline(meanE, '--'); yline(rmsE, ':'); yline(-rmsE, ':');

subplot(3,1,2);
plot(t, N, 'LineWidth', 1); grid on; ylabel('N [m]'); title('N Error');
yline(meanN, '--'); yline(rmsN, ':'); yline(-rmsN, ':');

subplot(3,1,3);
plot(t, U, 'LineWidth', 1); grid on; ylabel('U [m]'); xlabel('Epoch'); title('U Error');
yline(meanU, '--'); yline(rmsU, ':'); yline(-rmsU, ':');

%% 콘솔 요약 출력
fprintf('==== Horizontal & Vertical Error Stats ====\n');
fprintf('Mean   [E N U] = [%.3f  %.3f  %.3f] m\n', meanE, meanN, meanU);
fprintf('STD    [E N U] = [%.3f  %.3f  %.3f] m\n', stdE, stdN, stdU);
fprintf('RMS    [E N U] = [%.3f  %.3f  %.3f] m\n', rmsE, rmsN, rmsU);
fprintf('RMS2D  = %.3f m,   2DRMS = %.3f m\n', rms2D, d2rms);
fprintf('RMS3D  = %.3f m,   2DRMS = %.3f m\n', rms3D, d2rms);
fprintf('CEP50  = %.3f m,   CEP95 = %.3f m\n', cep50, cep95);




function LLA = parse_bestposa_lla(inputData)
%PARSE_BESTPOSA_LLA  Parse NovAtel #BESTPOSA lines to extract [lat lon hgt]
%   LLA = parse_bestposa_lla(filePath)            % from text file
%   LLA = parse_bestposa_lla(string_or_cellstr)   % from string / cellstr
%
% Output:
%   LLA: Nx3 double, columns = [lat_deg, lon_deg, hgt_m]
%
% Notes:
%   - 이 코드는 "#BESTPOSA" 라인 중 "SOL_COMPUTED,SINGLE,<lat>,<lon>,<hgt>" 패턴만 사용합니다.
%   - 불완전한 라인은 자동으로 건너뜁니다.

    % --- 1) 입력을 라인 셀 배열로 통일
    if ischar(inputData) || (isstring(inputData) && isscalar(inputData))
        % 파일 경로일 수도 있고, 단일 긴 문자열일 수도 있음
        if exist(char(inputData), 'file')
            % 파일에서 읽기
            txt = fileread(char(inputData));
            lines = regexp(txt, '\r\n|\n|\r', 'split');
        else
            % 단일 문자열을 여러 줄로 분할
            lines = regexp(char(inputData), '\r\n|\n|\r', 'split');
        end
    elseif isstring(inputData) || iscellstr(inputData)
        lines = cellstr(inputData);
    else
        error('Unsupported input type. Provide file path, string, or cellstr.');
    end

    % --- 2) 각 라인에서 lat, lon, hgt 추출
    LLA_list = [];  %#ok<NASGU>
    LLA_list = zeros(0,3);
    % 패턴: SOL_COMPUTED,SINGLE,<lat>,<lon>,<hgt>
    pat = 'SOL_COMPUTED,SINGLE,([^,]+),([^,]+),([^,]+)';

    for i = 1:numel(lines)
        line = strtrim(lines{i});
        if isempty(line) || ~contains(line, 'BESTPOSA')
            continue;
        end
        tok = regexp(line, pat, 'tokens', 'once');
        if ~isempty(tok)
            lat  = str2double(tok{1});
            lon  = str2double(tok{2});
            hgt  = str2double(tok{3});
            if all(~isnan([lat,lon,hgt]))
                LLA_list(end+1, :) = [lat, lon, hgt]; %#ok<AGROW>
            end
        end
    end

    LLA = LLA_list;
end
