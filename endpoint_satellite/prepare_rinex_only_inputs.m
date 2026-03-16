function [input_dir, nav_file, output_dir, selected_obs_file] = prepare_rinex_only_inputs(endpoint_root, varargin)
    data_dir = fullfile(endpoint_root, 'data');
    input_dir = fullfile(endpoint_root, 'test_input');
    output_dir = fullfile(endpoint_root, 'test_run_output');

    if ~exist(data_dir, 'dir')
        error('Data directory not found: %s', data_dir);
    end

    obs_override = '';
    nav_override = '';
    if nargin >= 2 && ~isempty(varargin{1})
        obs_override = varargin{1};
    end
    if nargin >= 3 && ~isempty(varargin{2})
        nav_override = varargin{2};
    end

    selected_obs_file = resolve_rinex_observation_file(data_dir, obs_override);
    nav_file = resolve_rinex_navigation_file(data_dir, nav_override);

    if exist(input_dir, 'dir')
        delete(fullfile(input_dir, '*'));
    else
        mkdir(input_dir);
    end

    [~, obs_name, ~] = fileparts(selected_obs_file);
    staged_obs_file = fullfile(input_dir, [obs_name '.25o']);
    copyfile(selected_obs_file, staged_obs_file, 'f');

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fprintf('RINEX-only observation file: %s\n', selected_obs_file);
    fprintf('RINEX-only navigation file: %s\n', nav_file);
    fprintf('RINEX-only staged observation file: %s\n', staged_obs_file);
end

function obs_file = resolve_rinex_observation_file(data_dir, obs_override)
    if ~isempty(obs_override)
        if exist(obs_override, 'file')
            obs_file = obs_override;
            return;
        end
        error('Observation file override not found: %s', obs_override);
    end

    candidates = find_rinex_candidate_files(data_dir);
    matches = {};
    for i = 1:numel(candidates)
        file_path = candidates{i};
        if is_rinex_observation_file(file_path)
            matches{end + 1} = file_path; %#ok<AGROW>
        end
    end

    if isempty(matches)
        error(['No observation RINEX file found under %s. ' ...
               'Place an observation file such as *_MO.rnx or *.25o in endpoint_satellite/data.'], data_dir);
    end

    obs_file = pick_most_recent_rinex(matches);
end

function nav_file = resolve_rinex_navigation_file(data_dir, nav_override)
    if ~isempty(nav_override)
        if exist(nav_override, 'file')
            nav_file = nav_override;
            return;
        end
        error('Navigation file override not found: %s', nav_override);
    end

    candidates = find_rinex_candidate_files(data_dir);
    matches = {};
    for i = 1:numel(candidates)
        file_path = candidates{i};
        if is_rinex_navigation_file(file_path)
            matches{end + 1} = file_path; %#ok<AGROW>
        end
    end

    if isempty(matches)
        error(['No navigation RINEX file found under %s. ' ...
               'Place a navigation file such as *_MN.rnx, *.25N, or *.25p in endpoint_satellite/data.'], data_dir);
    end

    nav_file = pick_most_recent_rinex(matches);
end

function file_list = find_rinex_candidate_files(root_dir)
    patterns = {'*.rnx', '*.RNX', '*.obs', '*.OBS', '*.25o', '*.25O', '*.25n', '*.25N', '*.25p', '*.25P', '*.nav', '*.NAV'};
    file_list = {};
    for i = 1:numel(patterns)
        items = dir(fullfile(root_dir, '**', patterns{i}));
        for j = 1:numel(items)
            if ~items(j).isdir
                file_list{end + 1} = fullfile(items(j).folder, items(j).name); %#ok<AGROW>
            end
        end
    end
    file_list = unique(file_list, 'stable');
end

function tf = is_rinex_observation_file(file_path)
    first_line = read_rinex_first_line(file_path);
    if contains(first_line, 'OBSERVATION DATA')
        tf = true;
        return;
    end

    name = lower(file_path);
    tf = contains(name, '_mo.') || endsWith(name, '.25o') || endsWith(name, '.obs');
end

function tf = is_rinex_navigation_file(file_path)
    first_line = read_rinex_first_line(file_path);
    if contains(first_line, 'NAVIGATION DATA') || contains(first_line, 'NAV DATA')
        tf = true;
        return;
    end

    name = lower(file_path);
    tf = contains(name, '_mn.') || contains(name, '_gn.') || endsWith(name, '.25n') || endsWith(name, '.25p') || endsWith(name, '.nav');
end

function first_line = read_rinex_first_line(file_path)
    fid = fopen(file_path, 'r');
    if fid < 0
        first_line = "";
        return;
    end

    cleanup = onCleanup(@() fclose(fid));
    first_line = string(fgetl(fid));
end

function selected_file = pick_most_recent_rinex(file_list)
    newest_time = -inf;
    newest_index = 1;

    for i = 1:numel(file_list)
        info = dir(file_list{i});
        current_time = info.datenum;
        if i == 1 || current_time > newest_time
            newest_time = current_time;
            newest_index = i;
        end
    end

    selected_file = file_list{newest_index};
end
