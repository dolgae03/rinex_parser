function summary_json = run_tsv_nav_endpoint_shared(tsv_file, nav_file, output_tsv)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    endpoint_root = fullfile(repo_root, 'endpoint_satellite');

    addpath(genpath(endpoint_root));
    addpath(genpath(repo_root));

    if nargin < 1 || isempty(tsv_file)
        error('A TSV input path is required.');
    end
    if nargin < 2 || isempty(nav_file)
        nav_file = '';
    end
    if nargin < 3 || isempty(output_tsv)
        output_tsv = '';
    end

    if isempty(output_tsv)
        [tsv_file, nav_file_resolved, output_tsv] = prepare_tsv_nav_inputs(endpoint_root, tsv_file, nav_file);
    else
        [tsv_file, nav_file_resolved, ~] = prepare_tsv_nav_inputs(endpoint_root, tsv_file, nav_file);
    end

    summary = fill_tsv_with_rinex_nav(tsv_file, nav_file_resolved, output_tsv);
    summary.input_tsv = tsv_file;
    summary.resolved_nav_file = nav_file_resolved;
    summary_json = jsonencode(summary);
end
