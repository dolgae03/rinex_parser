function summaries = process_and_export_rinex_files_endpoint(input_dir, nav_dir, output_dir, display_rows)
    if nargin < 4 || isempty(display_rows)
        display_rows = 10;
    end

    input_dir = normalize_text_input(input_dir, 'input_dir');
    nav_dir = normalize_text_input(nav_dir, 'nav_dir');
    output_dir = normalize_text_input(output_dir, 'output_dir');

    if ~exist(input_dir, 'dir')
        error('input_dir does not exist: %s', input_dir);
    end
    if ~exist(nav_dir, 'file')
        error('nav_dir does not exist: %s', nav_dir);
    end

    mat_output_dir = fullfile(output_dir, 'processed_mat');
    tsv_output_dir = fullfile(output_dir, 'satellite_tsv');

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    if ~exist(mat_output_dir, 'dir')
        mkdir(mat_output_dir);
    end
    if ~exist(tsv_output_dir, 'dir')
        mkdir(tsv_output_dir);
    end

    fprintf('=== Step 1/2: RINEX -> MAT ===\n');
    process_and_save_rinex_files(input_dir, nav_dir, mat_output_dir);

    fprintf('=== Step 2/2: MAT -> TSV with satellite positions ===\n');
    mat_files = find_all_files_recursive(mat_output_dir, '_processed.mat');
    summaries = [];

    for i = 1:numel(mat_files)
        mat_file = mat_files{i};
        relative_path = erase(mat_file, [mat_output_dir filesep]);
        [relative_folder, file_name, ~] = fileparts(relative_path);

        target_folder = fullfile(tsv_output_dir, relative_folder);
        if ~exist(target_folder, 'dir')
            mkdir(target_folder);
        end

        out_tsv = fullfile(target_folder, [file_name '.csv']);
        summary = export_processed_mat_to_tsv(mat_file, out_tsv, display_rows);
        if isempty(summaries)
            summaries = summary;
        else
            summaries(end + 1, 1) = summary; %#ok<AGROW>
        end

        fprintf('\n[%d/%d] %s\n', i, numel(mat_files), file_name);
        fprintf('  rows=%d, epochs=%d, output=%s\n', summary.row_count, summary.epoch_count, out_tsv);
        fprintf('  finite sv_pos rows=%d, finite sv_vel rows=%d, finite sv_clock rows=%d\n', ...
                summary.finite_sv_pos_rows, summary.finite_sv_vel_rows, summary.finite_sv_clock_rows);
        if ~isempty(summary.preview)
            disp(summary.preview);
        end

        if summary.finite_sv_pos_rows == 0
            error(['No finite satellite positions were produced for %s. ' ...
                   'This usually means the navigation file does not match the observation date ' ...
                   'or the selected constellation ephemerides are unavailable.'], file_name);
        end
    end
end

function text_value = normalize_text_input(value, name)
    if isstring(value) && isscalar(value)
        text_value = char(value);
        return;
    end

    if ischar(value)
        text_value = value;
        return;
    end

    error('%s must be a char array or string scalar.', name);
end

function file_list = find_all_files_recursive(root_dir, suffix)
    file_list = {};
    items = dir(root_dir);

    for i = 1:numel(items)
        item = items(i);
        current_path = fullfile(root_dir, item.name);

        if item.isdir
            if ~ismember(item.name, {'.', '..'})
                file_list = [file_list, find_all_files_recursive(current_path, suffix)]; %#ok<AGROW>
            end
        elseif endsWith(item.name, suffix)
            file_list{end + 1} = current_path; %#ok<AGROW>
        end
    end
end
