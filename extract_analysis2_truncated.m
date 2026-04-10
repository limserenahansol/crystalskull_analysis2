% EXTRACT: Baseline vs Injection with matched 10 min vs 10 min logic:
%   • Baseline trace: UNCHANGED (full ~10 min recording — no skip, no crop).
%   • Injection trace: LAST 10 minutes only (so long morphine/post sessions compare
%     to the same wall-clock length as baseline).
%
% Uses per-arm params: use_last_sec_inj only; baseline left full.
%
% Requires: extract_analysis2_core.m, extract_resolve_time_window_params.m,
%           extract_apply_time_window_T.m on the path.

clear; close all; clc;
tic;

%% --- Paths & parameters -------------------------------------------------
base_dir = '\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\CS1014-1b';
bl_path  = fullfile(base_dir, 'Baseline',  'mirror2path');
inj_path = fullfile(base_dir, 'Injection', 'mirror2path');
out_dir  = 'C:\Users\hsollim\Documents\MATLAB\MATLAB\2p\extract_figures_truncated';

inj_last_sec = 600;  % last 10 min of injection only (baseline stays full)

run_id = sprintf('CS1014-1b · mirror2path · BL_full_INJ_last%ds', inj_last_sec);

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
    'use_last_sec_inj', inj_last_sec);

extract_analysis2_core(bl_path, inj_path, out_dir, run_id, params);

fprintf('\nTotal wall time (script): %.1f s\n', toc);
