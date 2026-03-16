function nav_file = resolve_or_download_mixed_nav(tsv_file, local_nav_dir, cache_dir)
    validateattributes(tsv_file, {'char', 'string'}, {'scalartext'});
    validateattributes(local_nav_dir, {'char', 'string'}, {'scalartext'});
    validateattributes(cache_dir, {'char', 'string'}, {'scalartext'});

    tsv_file = char(tsv_file);
    local_nav_dir = char(local_nav_dir);
    cache_dir = char(cache_dir);

    day_infos = infer_nav_day_infos(tsv_file);
    requested_systems = infer_requested_systems(tsv_file);
    if ~exist(cache_dir, 'dir')
        mkdir(cache_dir);
    end

    daily_nav_files = cell(numel(day_infos), 1);
    for i = 1:numel(day_infos)
        daily_nav_files{i} = resolve_daily_nav_file(day_infos(i), local_nav_dir, cache_dir, requested_systems);
    end

    if numel(daily_nav_files) == 1
        nav_file = daily_nav_files{1};
    else
        nav_file = merge_daily_nav_files(daily_nav_files, cache_dir, day_infos);
    end
end

function requested_systems = infer_requested_systems(tsv_file)
    opts = detectImportOptions(tsv_file, 'FileType', 'text', 'Delimiter', '\t');
    opts.VariableNamingRule = 'preserve';
    if ~ismember('constellation', opts.VariableNames)
        requested_systems = 'GRECJ';
        return;
    end

    opts.SelectedVariableNames = {'constellation'};
    T = readtable(tsv_file, opts);
    values = unique(double(T.constellation(isfinite(double(T.constellation)))));

    requested_systems = '';
    for i = 1:numel(values)
        switch values(i)
            case 0
                requested_systems(end + 1) = 'G'; %#ok<AGROW>
            case 1
                requested_systems(end + 1) = 'E'; %#ok<AGROW>
            case 2
                requested_systems(end + 1) = 'C'; %#ok<AGROW>
            case 3
                requested_systems(end + 1) = 'R'; %#ok<AGROW>
            case 4
                requested_systems(end + 1) = 'J'; %#ok<AGROW>
            case 5
                error(['The TSV requests SBAS constellation data, but the current goGPS navigation parser ' ...
                       'fails on mixed-nav SBAS records. Wrapper-level filtering cannot safely enable SBAS.']);
            case 6
                requested_systems(end + 1) = 'I'; %#ok<AGROW>
        end
    end

    if isempty(requested_systems)
        requested_systems = 'GRECJ';
    end

    requested_systems = unique(requested_systems, 'stable');
end

function day_infos = infer_nav_day_infos(tsv_file)
    opts = detectImportOptions(tsv_file, 'FileType', 'text', 'Delimiter', '\t');
    opts.VariableNamingRule = 'preserve';

    selected_vars = intersect({'t_sec', 'gps_week', 'tow_sec'}, opts.VariableNames, 'stable');
    if isempty(selected_vars)
        error('The TSV must contain t_sec or gps_week/tow_sec to resolve a navigation date: %s', tsv_file);
    end
    opts.SelectedVariableNames = selected_vars;

    T = readtable(tsv_file, opts);
    times = extract_gps_datetimes(T);
    if isempty(times)
        error('No finite epoch time values were found in %s.', tsv_file);
    end

    day_starts = dateshift(times, 'start', 'day');
    unique_days = unique(day_starts, 'stable');

    day_infos = repmat(struct('start_time', [], 'year', 0, 'doy', 0, 'date_tag', ''), numel(unique_days), 1);
    for i = 1:numel(unique_days)
        day_infos(i).start_time = unique_days(i);
        day_infos(i).year = year(unique_days(i));
        day_infos(i).doy = day(unique_days(i), 'dayofyear');
        day_infos(i).date_tag = sprintf('%04d%03d', day_infos(i).year, day_infos(i).doy);
    end

    fprintf('Detected TSV time span: %s to %s (%d GPS day(s)).\n', ...
        char(min(times), 'yyyy-MM-dd HH:mm:ss'), ...
        char(max(times), 'yyyy-MM-dd HH:mm:ss'), ...
        numel(unique_days));
end

function times = extract_gps_datetimes(T)
    gps_epoch = datetime(1980, 1, 6, 0, 0, 0, 'TimeZone', 'UTC');

    if ismember('t_sec', T.Properties.VariableNames)
        t_sec = double(T.t_sec);
        valid = isfinite(t_sec);
        if any(valid)
            times = gps_epoch + seconds(t_sec(valid));
            return;
        end
    end

    if all(ismember({'gps_week', 'tow_sec'}, T.Properties.VariableNames))
        gps_week = double(T.gps_week);
        tow_sec = double(T.tow_sec);
        valid = isfinite(gps_week) & isfinite(tow_sec);
        if any(valid)
            times = gps_epoch + days(7 * gps_week(valid)) + seconds(tow_sec(valid));
            return;
        end
    end

    times = datetime.empty(0, 1);
end

function nav_file = resolve_daily_nav_file(day_info, local_nav_dir, cache_dir, requested_systems)
    search_roots = {local_nav_dir, cache_dir};
    candidate_paths = list_navigation_candidates(search_roots);

    for i = 1:numel(candidate_paths)
        [~, name, ext] = fileparts(candidate_paths{i});
        if contains(name, day_info.date_tag) && is_navigation_rinex(candidate_paths{i})
            nav_file = ensure_gogps_nav_filename(candidate_paths{i}, cache_dir, day_info, requested_systems);
            fprintf('Using local navigation file for %s: %s\n', day_info.date_tag, nav_file);
            return;
        end

        if strcmpi(ext, '.gz')
            decompressed = ensure_gunzip(candidate_paths{i}, cache_dir);
            [~, dec_name] = fileparts(decompressed);
            if contains(dec_name, day_info.date_tag) && is_navigation_rinex(decompressed)
                nav_file = ensure_gogps_nav_filename(decompressed, cache_dir, day_info, requested_systems);
                fprintf('Using local compressed navigation file for %s: %s\n', day_info.date_tag, nav_file);
                return;
            end
        end
    end

    nav_file = download_daily_nav_file(day_info, cache_dir, requested_systems);
end

function candidate_paths = list_navigation_candidates(search_roots)
    patterns = {'*.25N', '*.25n', '*.25P', '*.25p', '*_MN.rnx', '*_MN.RNX', '*.nav', '*.NAV', '*.gz', '*.GZ'};
    candidate_paths = {};
    for i = 1:numel(search_roots)
        if ~exist(search_roots{i}, 'dir')
            continue;
        end
        for j = 1:numel(patterns)
            items = dir(fullfile(search_roots{i}, '**', patterns{j}));
            for k = 1:numel(items)
                if ~items(k).isdir
                    candidate_paths{end + 1} = fullfile(items(k).folder, items(k).name); %#ok<AGROW>
                end
            end
        end
    end
    candidate_paths = unique(candidate_paths, 'stable');
end

function nav_file = download_daily_nav_file(day_info, cache_dir, requested_systems)
    base_url = sprintf('https://igs.bkg.bund.de/root_ftp/IGS/BRDC/%04d/%03d/', day_info.year, day_info.doy);
    file_names = { ...
        sprintf('BRDM00DLR_S_%s0000_01D_MN.rnx.gz', day_info.date_tag), ...
        sprintf('BRDC00WRD_S_%s0000_01D_MN.rnx.gz', day_info.date_tag), ...
        sprintf('BRDC00WRD_R_%s0000_01D_MN.rnx.gz', day_info.date_tag) ...
    };

    errors = cell(numel(file_names), 1);
    for i = 1:numel(file_names)
        gz_name = file_names{i};
        gz_path = fullfile(cache_dir, gz_name);
        uncompressed_path = fullfile(cache_dir, gz_name(1:end-3));
        if exist(uncompressed_path, 'file') && is_navigation_rinex(uncompressed_path)
            nav_file = ensure_gogps_nav_filename(uncompressed_path, cache_dir, day_info, requested_systems);
            fprintf('Using cached navigation file for %s: %s\n', day_info.date_tag, nav_file);
            return;
        end

        if exist(gz_path, 'file')
            try
                nav_file = ensure_gunzip(gz_path, cache_dir);
                if is_navigation_rinex(nav_file)
                    nav_file = ensure_gogps_nav_filename(nav_file, cache_dir, day_info, requested_systems);
                    fprintf('Using cached compressed navigation file for %s: %s\n', day_info.date_tag, nav_file);
                    return;
                end
            catch me
                errors{i} = sprintf('%s (cached gunzip failed: %s)', gz_name, me.message);
            end
        end

        url = [base_url gz_name];
        fprintf('Downloading mixed navigation file: %s\n', url);
        try
            websave(gz_path, url, weboptions('Timeout', 60));
            nav_file = ensure_gunzip(gz_path, cache_dir);
            if ~is_navigation_rinex(nav_file)
                error('Downloaded file is not a navigation RINEX: %s', nav_file);
            end
            nav_file = ensure_gogps_nav_filename(nav_file, cache_dir, day_info, requested_systems);
            return;
        catch me
            errors{i} = sprintf('%s -> %s', url, me.message);
            if exist(gz_path, 'file')
                delete(gz_path);
            end
        end
    end

    error('Unable to resolve a mixed navigation file for %s. Tried:%s\n- %s', ...
        day_info.date_tag, newline, strjoin(errors, [newline '- ']));
end

function nav_file = ensure_gunzip(gz_path, target_dir)
    if ~exist(gz_path, 'file')
        error('Compressed navigation file not found: %s', gz_path);
    end

    nav_file = gz_path;
    if ~endsWith(lower(gz_path), '.gz')
        return;
    end

    output_files = gunzip(gz_path, target_dir);
    if isempty(output_files)
        error('gunzip produced no output for %s.', gz_path);
    end
    nav_file = output_files{1};
end

function normalized_file = ensure_gogps_nav_filename(nav_file, cache_dir, day_info, requested_systems)
    if ~is_mixed_navigation_rinex(nav_file)
        normalized_file = nav_file;
        return;
    end

    two_digit_year = mod(day_info.year, 100);
    canonical_name = sprintf('brdm%03d0_%s.%02dp', day_info.doy, lower(requested_systems), two_digit_year);
    normalized_file = fullfile(cache_dir, canonical_name);

    if exist(normalized_file, 'file')
        return;
    end

    write_filtered_mixed_nav(nav_file, normalized_file, requested_systems);
end

function merged_file = merge_daily_nav_files(daily_nav_files, cache_dir, day_infos)
    merged_file = fullfile(cache_dir, ...
        sprintf('BRDC_MERGED_%s_%s_01D_MN.25p', day_infos(1).date_tag, day_infos(end).date_tag));
    if exist(merged_file, 'file')
        return;
    end

    out_fid = fopen(merged_file, 'w');
    if out_fid < 0
        error('Unable to create merged navigation file: %s', merged_file);
    end
    cleanup = onCleanup(@() fclose(out_fid));

    for i = 1:numel(daily_nav_files)
        in_fid = fopen(daily_nav_files{i}, 'r');
        if in_fid < 0
            error('Unable to open navigation file for merge: %s', daily_nav_files{i});
        end

        write_line = (i == 1);
        while ~feof(in_fid)
            line = fgetl(in_fid);
            if ~ischar(line)
                break;
            end

            if ~write_line
                if contains(line, 'END OF HEADER')
                    write_line = true;
                end
                continue;
            end

            fprintf(out_fid, '%s\n', line);
        end
        fclose(in_fid);
    end
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

function tf = is_mixed_navigation_rinex(file_path)
    fid = fopen(file_path, 'r');
    if fid < 0
        tf = false;
        return;
    end

    cleanup = onCleanup(@() fclose(fid));
    first_line = char(fgetl(fid));
    tf = contains(first_line, 'NAVIGATION DATA') && contains(first_line, 'M');
end

function write_filtered_mixed_nav(source_file, target_file, requested_systems)
    in_fid = fopen(source_file, 'rt');
    if in_fid < 0
        error('Unable to open source navigation file: %s', source_file);
    end

    out_fid = fopen(target_file, 'wt');
    if out_fid < 0
        fclose(in_fid);
        error('Unable to create filtered navigation file: %s', target_file);
    end

    in_cleanup = onCleanup(@() fclose(in_fid));
    out_cleanup = onCleanup(@() fclose(out_fid));

    while ~feof(in_fid)
        line = fgetl(in_fid);
        if ~ischar(line)
            break;
        end
        fprintf(out_fid, '%s\n', line);
        if contains(line, 'END OF HEADER')
            break;
        end
    end

    while ~feof(in_fid)
        line1 = fgetl(in_fid);
        if ~ischar(line1)
            break;
        end
        if isempty(strtrim(line1))
            continue;
        end

        sys_char = line1(1);
        n_lines = get_nav_record_line_count(sys_char);
        record_lines = cell(n_lines, 1);
        record_lines{1} = pad_nav_line(line1);

        for i = 2:n_lines
            next_line = fgetl(in_fid);
            if ~ischar(next_line)
                error('Unexpected EOF while filtering navigation record starting with %s', line1);
            end
            record_lines{i} = pad_nav_line(next_line);
        end

        if contains(requested_systems, sys_char)
            for i = 1:n_lines
                fprintf(out_fid, '%s\n', record_lines{i});
            end
        end
    end
end

function n_lines = get_nav_record_line_count(sys_char)
    switch sys_char
        case {'R', 'S'}
            n_lines = 4;
        otherwise
            n_lines = 8;
    end
end

function line = pad_nav_line(line)
    if numel(line) < 80
        line = [line repmat(' ', 1, 80 - numel(line))];
    end
end
