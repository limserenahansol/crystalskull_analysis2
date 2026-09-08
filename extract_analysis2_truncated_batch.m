function extract_analysis2_truncated_batch()
%EXTRACT_ANALYSIS2_TRUNCATED_BATCH  Same as extract_analysis2_batch, but:
%   Injection = last N seconds where N = baseline duration after any BL time crop (per mouse).
%   Default: full baseline file + matched injection tail (900 s BL → last 900 s INJ; 800 → last 800 s).
%   Optional: set use_first_sec_bl so N = that cropped length (e.g. first 600 s BL → last 600 s INJ).
%
%   extract_analysis2_truncated_batch
%
% Requires: extract_analysis2_core.m, extract_resolve_time_window_params.m,
%           extract_apply_time_window_T.m, extract_resolve_extract_files.m

tic;

%% ========== CONFIG ======================================================
base_dir = '\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\meltinonly';

mirror_subfolders = {'mirror1path', 'mirror2path'};
baseline_folders = {'Baseline'};

out_root = 'C:\Users\hsollim\Documents\MATLAB\MATLAB\2p\extract_figures_batch_meltinonly_truncated';

mouse_ids_only = {};

% Optional: uncomment to use only the first N seconds of baseline; injection tail length N matches that.
% bl_cap_sec = 600;
% use_first_sec_bl in params below = bl_cap_sec

params = struct( ...
    'z_thresh', 2, ...
    'n_boot', 1000, ...
    'time_bin', 1, ...
    'n_pc', 10, ...
    'n_corr_sample', 500, ...
    'frame_px', 2304, ...
    'max_show', 2000, ...
    'skip_first_sec_bl', 0, ...
    'use_last_sec_bl', inf, ...
    'skip_first_sec_inj', 0, ...
    'use_last_sec_inj', inf, ...
    'match_inj_last_to_baseline_dur', true);
%% ========================================================================

close all;
clc;

if isempty(mouse_ids_only)
    mouse_ids_only = discover_mouse_folders(base_dir);
end

if isempty(mouse_ids_only)
    error('No mouse folders found under base_dir with Baseline and Injection: %s', base_dir);
end

if exist(out_root, 'dir') ~= 7
    mkdir(out_root);
end

n_ok = 0;
n_skip = 0;
n_fail = 0;

for mi = 1:numel(mouse_ids_only)
    mouse_id = mouse_ids_only{mi};
    mouse_root = fullfile(base_dir, mouse_id);

    for bi = 1:numel(baseline_folders)
        bl_folder = baseline_folders{bi};

        for ki = 1:numel(mirror_subfolders)
            mirror_name = mirror_subfolders{ki};
            bl_path  = fullfile(mouse_root, bl_folder, mirror_name);
            inj_path = fullfile(mouse_root, 'Injection', mirror_name);

            [ready, skip_reason] = session_pair_ready(bl_path, inj_path);
            if ~ready
                % Log line names the session key; failure may be on INJ even when BL folder is complete.
                fprintf('SKIP  %s / %s / %s\n  Reason: %s\n  BL:  %s\n  INJ: %s\n', ...
                    mouse_id, bl_folder, mirror_name, skip_reason, bl_path, inj_path);
                n_skip = n_skip + 1;
                continue;
            end

            out_tag = sprintf('%s_%s_%s_INJlastMatchBLdur', mouse_id, bl_folder, mirror_name);
            out_dir = fullfile(out_root, out_tag);
            run_id  = sprintf('%s · %s · %s · INJ last window = baseline duration', mouse_id, bl_folder, mirror_name);

            fprintf('\n-------- RUN %s --------\n', out_tag);
            try
                extract_analysis2_core(bl_path, inj_path, out_dir, run_id, params);
                n_ok = n_ok + 1;
            catch ME
                warning('extract_analysis2_truncated_batch:RunFailed', 'FAILED %s: %s', out_tag, ME.message);
                n_fail = n_fail + 1;
            end
            close all;
        end
    end
end

fprintf('\n======== BATCH DONE (truncated) ========\n');
fprintf('  OK:    %d\n  Skip:  %d\n  Fail:  %d\n', n_ok, n_skip, n_fail);
fprintf('  Total wall time: %.1f s (%.1f min)\n', toc, toc/60);

end

function mouse_ids = discover_mouse_folders(base_dir)
mouse_ids = {};
if exist(base_dir, 'dir') ~= 7
    return;
end
d = dir(base_dir);
for i = 1:numel(d)
    if ~d(i).isdir || strcmp(d(i).name, '.') || strcmp(d(i).name, '..')
        continue;
    end
    mp = fullfile(base_dir, d(i).name);
    if exist(fullfile(mp, 'Baseline'), 'dir') == 7 && exist(fullfile(mp, 'Injection'), 'dir') == 7
        mouse_ids{end+1} = d(i).name; %#ok<AGROW>
    end
end
end

function [tf, reason] = session_pair_ready(bl_path, inj_path)
%SESSION_PAIR_READY  True only if BL+INJ mirror folders exist, cell_map.png present, mats resolved.
tf = false;
reason = '';
if exist(bl_path, 'dir') ~= 7
    reason = sprintf('Baseline mirror folder missing: %s', bl_path);
    return;
end
if exist(inj_path, 'dir') ~= 7
    reason = sprintf('Injection mirror folder missing: %s', inj_path);
    return;
end
if exist(fullfile(bl_path, 'cell_map.png'), 'file') ~= 2
    reason = sprintf('Baseline mirror has no cell_map.png (%s)', bl_path);
    return;
end
if exist(fullfile(inj_path, 'cell_map.png'), 'file') ~= 2
    reason = sprintf('Injection mirror has no cell_map.png (%s)', inj_path);
    return;
end
Rb = extract_resolve_extract_files(bl_path);
Ri = extract_resolve_extract_files(inj_path);
miss_b = missing_extract_labels(Rb);
miss_i = missing_extract_labels(Ri);
if ~isempty(miss_b)
    reason = sprintf('Baseline mats unresolved (%s): need %s', bl_path, miss_b);
    return;
end
if ~isempty(miss_i)
    reason = sprintf('Injection mats unresolved (%s): need %s', inj_path, miss_i);
    return;
end
tf = true;
reason = '';
end

function s = missing_extract_labels(R)
% Names expected by extract_resolve_extract_files (frr, ds_ext, summary).
need = {'frr', 'ds_ext', 'summary'};
miss = {};
for k = 1:numel(need)
    if isempty(R.(need{k}))
        miss{end+1} = need{k}; %#ok<AGROW>
    end
end
if isempty(miss)
    s = '';
else
    s = strjoin(miss, ', ');
end
end
