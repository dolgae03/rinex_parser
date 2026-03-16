repo_root = fileparts(fileparts(mfilename('fullpath')));
endpoint_root = fullfile(repo_root, 'endpoint_satellite');

addpath(genpath(endpoint_root));
addpath(genpath(repo_root));

tsv_file = getenv('SAT_ENDPOINT_TSV_FILE');
nav_file = getenv('SAT_ENDPOINT_NAV_FILE');
output_tsv = getenv('SAT_ENDPOINT_OUTPUT_TSV');

if isempty(tsv_file) || isempty(nav_file) || isempty(output_tsv)
    [tsv_file, nav_file, output_tsv] = prepare_tsv_nav_inputs(endpoint_root);
end

summary = fill_tsv_with_rinex_nav(tsv_file, nav_file, output_tsv); %#ok<NASGU>
disp(summary);
