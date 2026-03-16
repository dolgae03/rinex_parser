repo_root = 'C:\Users\mskim\Desktop\workspace\goGPS_loadRinex';
endpoint_root = fullfile(repo_root, 'endpoint_satellite');

addpath(genpath(endpoint_root));
addpath(genpath(repo_root));

[input_dir, nav_file, output_dir, selected_obs_file] = prepare_rinex_only_inputs(endpoint_root);
display_rows = 5;

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

summaries = process_and_export_rinex_files_endpoint(input_dir, nav_file, output_dir, display_rows); %#ok<NASGU>
disp(['run_rinex_only_current_test observation: ' selected_obs_file]);
disp('run_rinex_only_current_test finished');
