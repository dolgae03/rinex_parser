
convert_mat_182_fields_to_csv("./data/correction_20230523_1916_pixel5.mat","./data/csv/20230523_1916_pixel5","t1s_UTC"); % 시간 필드 지정
convert_mat_182_fields_to_csv("./data/correction_20230519_2010_pixel5.mat","./data/csv/20230519_2010_pixel5","t1s_UTC"); % 시간 필드 지정
convert_mat_182_fields_to_csv("./data/correction_20230509_2132_pixel5.mat","./data/csv/20230509_2132_pixel5","t1s_UTC"); % 시간 필드 지정
convert_mat_182_fields_to_csv("./data/correction_20220804_2007_pixel5.mat","./data/csv/20220804_2007_pixel5","t1s_UTC"); % 시간 필드 지정
convert_mat_182_fields_to_csv("./data/correction_20220126_2002_pixel5.mat","./data/csv/20220126_2002_pixel5","t1s_UTC"); % 시간 필드 지정

function convert_mat_182_fields_to_csv(mat_filename, out_dir, time_field)
% convert_mat_182_fields_to_csv
%   MAT 파일의 모든 필드 중에서 "행 또는 열 차원 중 하나가 182"인 numeric 데이터들을 찾아
%   각 필드명으로 CSV를 생성합니다.
%   CSV 앞부분에는 UTC_Time(ISO8601, UTC), UnixMillis(UTC epoch ms)를 공통 컬럼으로 추가합니다.
%
% 사용 예:
%   convert_mat_182_fields_to_csv("matlab.mat","./out");                   % 자동 시간 추정
%   convert_mat_182_fields_to_csv("matlab.mat","./out","t1s_UTC");         % 시간 필드 지정
%
% 인자:
%   mat_filename : 입력 mat 경로
%   out_dir      : CSV 저장 폴더
%   time_field   : (옵션) 시간 벡터 필드명. 생략/빈 문자열이면 자동 탐지

    if nargin < 2 || isempty(out_dir), out_dir = "."; end
    if ~isfolder(out_dir), mkdir(out_dir); end

    S = load(mat_filename);
    fns = fieldnames(S);

    % --- 시간 벡터 준비 ---
    [posix_sec, time_src_name] = get_time_vector_posix(S, nargin>=3 && ~isempty(time_field), time_field);

    % 공통 시간 컬럼
    t_millis = round(posix_sec * 1000);           % double로 저장(엑셀 호환)
    utc_time = datetime(posix_sec,'ConvertFrom','posixtime','TimeZone','UTC');

    fprintf('Time source: %s (N=%d)\n', time_src_name, numel(posix_sec));

    % --- 모든 182-채널 필드 탐색 & 저장 ---
    saved = 0;
    for i = 1:numel(fns)
        name = fns{i};
        if strcmp(name, time_src_name), continue; end     % 시간 필드는 스킵
        val = S.(name);

        if ~isnumeric(val), continue; end
        sz = size(val);

        % 2D만 처리 (3D 이상은 스킵)
        if numel(sz) ~= 2, continue; end

        % 182 차원을 포함하는가?
        has182 = any(sz == 182);
        if ~has182, continue; end

        % 어느 축이 time과 맞는지 판단
        N = numel(posix_sec);
        [r,c] = size(val);

        % 목표 형태: (N_time x 182)
        orient_ok = false;
        if r == N && c == 182
            C = val;
            orient_ok = true;
        elseif r == 182 && c == N
            C = val.';   % 전치
            orient_ok = true;
        elseif c == 182 && r == 1 % 단일 샘플(1x182) → 시간 1개로 저장 허용
            C = val;
            if N ~= 1
                warning('필드 %s: 데이터는 1x182인데 시간 길이 N=%d와 불일치 → 스킵', name, N);
                continue;
            end
            orient_ok = true;
        elseif r == 182 && c == 1 % 단일 샘플(182x1)
            C = val.';  % 1x182로
            if N ~= 1
                warning('필드 %s: 데이터는 182x1인데 시간 길이 N=%d와 불일치 → 스킵', name, N);
                continue;
            end
            orient_ok = true;
        else
            % 시간 길이가 맞지 않으면 스킵
            if r == N
                warning('필드 %s: 열 개수 %d가 182가 아님 → 스킵', name, c);
            elseif c == N
                warning('필드 %s: 행 개수 %d가 182가 아님 → 스킵', name, r);
            else
                warning('필드 %s: 시간 길이(%d)와 어느 축도 일치하지 않음(%dx%d) → 스킵', name, N, r, c);
            end
            continue;
        end

        if ~orient_ok, continue; end

        % 변수명 corr_001 ~ corr_182
        nCols = size(C,2);
        if nCols ~= 182
            warning('필드 %s: 182 채널이 아님(현재 %d). 그래도 저장합니다.', name, nCols);
        end
        corrVarNames = arrayfun(@(k) sprintf('corr_%03d', k), 1:nCols, 'UniformOutput', false);

        T = table(utc_time, t_millis, 'VariableNames', {'UTC_Time','UnixMillis'});
        T = [T, array2table(double(C), 'VariableNames', corrVarNames)];

        csv_path = fullfile(out_dir, sprintf('%s.csv', name));
        writetable(T, csv_path);

        fprintf('✅ [%s] 저장: %s (행:%d, 채널:%d)\n', name, csv_path, size(C,1), size(C,2));
        saved = saved + 1;
    end

    if saved == 0
        fprintf('⚠️ 182-채널 형태를 가진 필드를 찾지 못했습니다.\n');
    else
        fprintf('🎉 총 %d개 필드를 CSV로 저장 완료.\n', saved);
    end
end

% ===== 내부 유틸 =====
function [posix_sec, src_name] = get_time_vector_posix(S, has_user_choice, time_field)
    % 1) 사용자가 명시한 시간 필드 우선
    if has_user_choice
        assert(isfield(S, time_field), '지정한 시간 필드(%s)가 MAT에 없습니다.', time_field);
        posix_sec = normalize_to_posix_seconds(S.(time_field));
        src_name = time_field;
        return;
    end

    % 2) 흔한 후보 우선순위
    candidates = {'t1s_UTC','ti','utc_time','UTC_Time','time','t','posix','unix','UnixMillis'};
    for k = 1:numel(candidates)
        nm = candidates{k};
        if isfield(S, nm)
            try
                posix_sec = normalize_to_posix_seconds(S.(nm));
                src_name = nm;
                return;
            catch
                % 실패 시 다음 후보
            end
        end
    end

    % 3) 자동 탐색: 단조 증가하는 벡터형(datetime/numeric) 찾기
    fns = fieldnames(S);
    best = [];
    bestname = '';
    for i = 1:numel(fns)
        v = S.(fns{i});
        try
            if isvector(v) && (isnumeric(v) || isdatetime(v))
                p = normalize_to_posix_seconds(v);
                if numel(p) >= 1 && is_monotonic_increasing(p)
                    best = p; bestname = fns{i};
                    break;
                end
            end
        catch
            % skip
        end
    end
    if ~isempty(best)
        posix_sec = best;
        src_name = bestname;
        return;
    end

    error('시간 벡터를 찾을 수 없습니다. time_field를 지정하세요 (예: "t1s_UTC").');
end

function tf = is_monotonic_increasing(x)
    x = x(:);
    dx = diff(x);
    tf = all(isfinite(x)) && all(dx >= 0) && any(dx > 0);
end

function posix_sec = normalize_to_posix_seconds(ti)
% 다양한 ti 입력을 안전하게 POSIX seconds(UTC)로 변환
% - datetime(UTC/naive) → posixtime
% - datenum(double)     → datetime(...,'ConvertFrom','datenum') → posixtime
% - numeric POSIX seconds / milliseconds / microseconds / nanoseconds 자동 판별

    if isdatetime(ti)
        if isempty(ti.TimeZone)
            ti.TimeZone = 'UTC';
        else
            ti = datetime(ti,'TimeZone','UTC');
        end
        posix_sec = posixtime(ti);
        posix_sec = double(posix_sec(:));

    elseif isnumeric(ti)
        ti = double(ti(:));
        mx = max(ti);

        % 스케일 판별:
        %   POSIX s (2020~2035): ~1.6e9~2.1e9
        %   ms: ~1.6e12
        %   us: ~1.6e15
        %   ns: ~1.6e18
        if mx > 1e15       % ns
            posix_sec = ti / 1e9;
        elseif mx > 1e12   % us
            posix_sec = ti / 1e6;
        elseif mx > 1e10   % ms
            posix_sec = ti / 1e3;
        elseif mx > 7e5    % datenum (일)
            dt = datetime(ti,'ConvertFrom','datenum','TimeZone','UTC');
            posix_sec = posixtime(dt);
            posix_sec = double(posix_sec(:));
        else               % s
            posix_sec = ti;
        end

    else
        error('지원하지 않는 시간 타입: %s', class(ti));
    end
end
