repo_root = fileparts(fileparts(mfilename('fullpath')));
endpoint_root = fullfile(repo_root, 'endpoint_satellite');

addpath(genpath(endpoint_root));
addpath(genpath(repo_root));

input_dir = getenv('SAT_ENDPOINT_INPUT_DIR');
nav_file = getenv('SAT_ENDPOINT_NAV_FILE');
output_dir = getenv('SAT_ENDPOINT_OUTPUT_DIR');
display_rows_text = getenv('SAT_ENDPOINT_DISPLAY_ROWS');

if isempty(input_dir) || isempty(nav_file) || isempty(output_dir)
    [input_dir, nav_file, output_dir] = prepare_rinex_only_inputs(endpoint_root);
end

if isempty(display_rows_text)
    display_rows = 10;
else
    display_rows = str2double(display_rows_text);
    if isnan(display_rows) || display_rows < 1
        display_rows = 10;
    end
end

process_and_export_rinex_files_endpoint(input_dir, nav_file, output_dir, display_rows);
