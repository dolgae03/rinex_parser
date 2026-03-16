% 데이터 파일 경로
data1_path = "D:\과제\삼성 과제\Data\Processed\1차년도 Data 수집\BMHR21030714D_2025-05-07_17-21-04_processed.mat";
data2_path = "C:\Users\mskim\Desktop\workspace\samsung_analyze_tool\data\obs_real\29742_gnss_log_2025_05_07_17_22_10.mat";

% 데이터 로드
data1 = load(data1_path);
data2 = load(data2_path);

% time_sync 함수 호출
[sync_data1, sync_data2] = timesync(data1, data2);

% 저장 경로 생성 함수
function aligned_path = make_aligned_path(original_path)
    [folder, name, ~] = fileparts(original_path);
    aligned_path = fullfile(folder, name + "_aligned.mat");
end

% 동기화된 데이터 저장 경로 생성
output_path1 = make_aligned_path(data1_path);
output_path2 = make_aligned_path(data2_path);

% 저장
save(output_path1, '-struct', 'sync_data1');
save(output_path2, '-struct', 'sync_data2');

disp('동기화된 데이터가 성공적으로 저장되었습니다.');
