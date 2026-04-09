function extract_analysis2_batch()
%EXTRACT_ANALYSIS2_BATCH  Batch Baseline vs Injection for all mice × mirrors.
%
%   extract_analysis2_batch
%
% Layout:
%   base_dir/
%     CS1014-1b/
%       Baseline/mirror1path/   (M_moco_frr.mat, M_moco_ds_ext.mat, M_summary.mat, cell_map.png)
%       Injection/mirror1path/
%       Baseline/mirror2path/
%       Injection/mirror2path/
%
% Outputs: fullfile(out_root, '<mouseID>_<mirrorName>')/
%
% Edit CONFIG below, then run: extract_analysis2_batch
%
% Requires: extract_analysis2_core.m on the MATLAB path.

tic;

%% ========== CONFIG ======================================================
base_dir = '\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon';

mirror_subfolders = {'mirror1path', 'mirror2path'};

out_root = 'C:\Users\hsollim\Documents\MATLAB\MATLAB\2p\extract_figures_batch';

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

req_files = {'M_moco_frr.mat', 'M_moco_ds_ext.mat', 'M_summary.mat', 'cell_map.png'};

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

    for ki = 1:numel(mirror_subfolders)
        mirror_name = mirror_subfolders{ki};
        bl_path  = fullfile(mouse_root, 'Baseline',  mirror_name);
        inj_path = fullfile(mouse_root, 'Injection', mirror_name);

        if ~data_ready(bl_path, inj_path, req_files)
            fprintf('SKIP  %s / %s  (missing files under BL or INJ)\n', mouse_id, mirror_name);
            n_skip = n_skip + 1;
            continue;
        end

        out_tag = sprintf('%s_%s', mouse_id, mirror_name);
        out_dir = fullfile(out_root, out_tag);
        run_id  = sprintf('%s · %s', mouse_id, mirror_name);

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

function tf = data_ready(bl_path, inj_path, req_files)
tf = false;
if exist(bl_path, 'dir') ~= 7 || exist(inj_path, 'dir') ~= 7
    return;
end
for j = 1:numel(req_files)
    f = req_files{j};
    if exist(fullfile(bl_path, f), 'file') ~= 2 || exist(fullfile(inj_path, f), 'file') ~= 2
        return;
    end
end
tf = true;
end
