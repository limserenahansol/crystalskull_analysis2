function extract_analysis2_batch()
%EXTRACT_ANALYSIS2_BATCH  Batch Baseline vs Injection for all mice × mirrors.
%
%   extract_analysis2_batch
%
% Layout:
%   base_dir/
%     1014m01b/
%       Baseline/mirror1path/   (*_moco_frr.mat, *_moco_ds_pre_ext.mat, *_summary.mat, cell_map.png)
%       Injection/mirror1path/
%     1014m02w/
%       Baseline/ ...  and optionally Baseline2/ ...
%       Injection/ ...
%
% Outputs: fullfile(out_root, '<mouseID>_<baselineFolder>_<mirrorName>')/
%
% Edit CONFIG below, then run: extract_analysis2_batch
%
% Requires: extract_analysis2_core.m on the MATLAB path.

tic;

%% ========== CONFIG ======================================================
base_dir = '\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\morphine_20260408';

mirror_subfolders = {'mirror1path', 'mirror2path'};

% Condition folder names on the baseline side (each paired with Injection/<mirror>)
baseline_folders = {'Baseline'};
% e.g. 1014m02w with a second baseline session:
% baseline_folders = {'Baseline', 'Baseline2'};

out_root = 'C:\Users\hsollim\Documents\MATLAB\MATLAB\2p\extract_figures_batch_morphin';

% If empty: auto-discover subfolders of base_dir that contain Baseline + Injection.
mouse_ids_only = {};

params = struct( ...
    'z_thresh', 2, ...
    'n_boot', 1000, ...
    'time_bin', 1, ...
    'n_pc', 10, ...
    'n_corr_sample', 500, ...
    'frame_px', 2304, ...
    'max_show', 2000);
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

            if ~session_pair_ready(bl_path, inj_path)
                fprintf('SKIP  %s / %s / %s  (missing folder, cell_map.png, or resolved .mat files)\n', ...
                    mouse_id, bl_folder, mirror_name);
                n_skip = n_skip + 1;
                continue;
            end

            out_tag = sprintf('%s_%s_%s', mouse_id, bl_folder, mirror_name);
            out_dir = fullfile(out_root, out_tag);
            run_id  = sprintf('%s · %s · %s', mouse_id, bl_folder, mirror_name);

            fprintf('\n-------- RUN %s --------\n', out_tag);
            try
                extract_analysis2_core(bl_path, inj_path, out_dir, run_id, params);
                n_ok = n_ok + 1;
            catch ME
                warning('extract_analysis2_batch:RunFailed', 'FAILED %s: %s', out_tag, ME.message);
                n_fail = n_fail + 1;
            end
            close all;
        end
    end
end

fprintf('\n======== BATCH DONE ========\n');
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

function tf = session_pair_ready(bl_path, inj_path)
tf = false;
if exist(bl_path, 'dir') ~= 7 || exist(inj_path, 'dir') ~= 7
    return;
end
if exist(fullfile(bl_path, 'cell_map.png'), 'file') ~= 2
    return;
end
if exist(fullfile(inj_path, 'cell_map.png'), 'file') ~= 2
    return;
end
Rb = extract_resolve_extract_files(bl_path);
Ri = extract_resolve_extract_files(inj_path);
tf = ~isempty(Rb.frr) && ~isempty(Rb.ds_ext) && ~isempty(Rb.summary) ...
    && ~isempty(Ri.frr) && ~isempty(Ri.ds_ext) && ~isempty(Ri.summary);
end
