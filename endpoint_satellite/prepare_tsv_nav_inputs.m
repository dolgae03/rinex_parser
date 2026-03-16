function [tsv_file, nav_file, output_tsv] = prepare_tsv_nav_inputs(endpoint_root, varargin)
    data_dir = fullfile(endpoint_root, 'data');
    output_dir = fullfile(endpoint_root, 'tsv_nav_output');

    if ~exist(data_dir, 'dir')
        error('Data directory not found: %s', data_dir);
    end

    tsv_override = '';
    nav_override = '';
    if nargin >= 2 && ~isempty(varargin{1})
        tsv_override = varargin{1};
    end
    if nargin >= 3 && ~isempty(varargin{2})
        nav_override = varargin{2};
    end

    tsv_file = resolve_tsv_file(data_dir, tsv_override);
    nav_file = resolve_nav_file(data_dir, nav_override);

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    [~, name, ext] = fileparts(tsv_file);
    output_tsv = fullfile(output_dir, [name '_filled' ext]);

    fprintf('TSV input file: %s\n', tsv_file);
    fprintf('Navigation RINEX file: %s\n', nav_file);
    fprintf('Output TSV file: %s\n', output_tsv);
end

function tsv_file = resolve_tsv_file(data_dir, tsv_override)
    if ~isempty(tsv_override)
        if exist(tsv_override, 'file')
            tsv_file = tsv_override;
            return;
        end
        error('TSV override not found: %s', tsv_override);
    end

    matches = dir(fullfile(data_dir, '**', '*.tsv'));
    if isempty(matches)
        error('No TSV input found under %s.', data_dir);
    end

    tsv_file = fullfile(matches(1).folder, matches(1).name);
end

function nav_file = resolve_nav_file(data_dir, nav_override)
    if ~isempty(nav_override)
        if exist(nav_override, 'file')
            nav_file = nav_override;
            return;
        end
        error('Navigation override not found: %s', nav_override);
    end

    patterns = {'*.25N', '*.25n', '*.25P', '*.25p', '*_MN.rnx', '*_MN.RNX', '*.nav', '*.NAV'};
    candidates = {};
    for i = 1:numel(patterns)
        items = dir(fullfile(data_dir, '**', patterns{i}));
        for j = 1:numel(items)
            if ~items(j).isdir
                file_path = fullfile(items(j).folder, items(j).name);
                if is_navigation_rinex(file_path)
                    candidates{end + 1} = file_path; %#ok<AGROW>
                end
            end
        end
    end

    if isempty(candidates)
        error('No navigation RINEX file found under %s.', data_dir);
    end

    nav_file = candidates{1};
end

function tf = is_navigation_rinex(file_path)
    fid = fopen(file_path, 'r');
    if fid < 0
        tf = false;
        return;
    end

    cleanup = onCleanup(@() fclose(fid));
    first_line = string(fgetl(fid));
    tf = contains(first_line, 'NAVIGATION DATA') || contains(first_line, 'NAV DATA');
end
