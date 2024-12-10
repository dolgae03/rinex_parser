% 데이터 파일 경로
data1_path = "D:\과제\삼성 과제\Data\Processed\baseline.mat";
data2_path = "D:\과제\삼성 과제\Data\Processed\rover.mat";

% 데이터 로드
data1 = load(data1_path);
data2 = load(data2_path);

% time_sync 함수 호출
[sync_data1, sync_data2] = timesync(data1, data2);

% 동기화된 데이터 저장 경로
output_path1 = "D:\과제\삼성 과제\Data\Processed\sync_baseline.mat";
output_path2 = "D:\과제\삼성 과제\Data\Processed\sync_rover.mat";

% 저장
save(output_path1, '-struct', 'sync_data1');
save(output_path2, '-struct', 'sync_data2');

disp('동기화된 데이터가 성공적으로 저장되었습니다.');
