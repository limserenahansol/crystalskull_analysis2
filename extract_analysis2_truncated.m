% EXTRACT Results: Baseline vs Injection — same as extract_analysis2.m, but after
% loading traces we DROP the first 10 minutes, then keep ONLY the last 10 minutes
% of the remainder (per condition, using each side's own frame rate).
%
% Requires: extract_analysis2_core.m, extract_apply_time_window_T.m on the path.

clear; close all; clc;
tic;

%% --- Paths & parameters -------------------------------------------------
base_dir = '\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\CS1014-1b';
bl_path  = fullfile(base_dir, 'Baseline',  'mirror2path');
inj_path = fullfile(base_dir, 'Injection', 'mirror2path');
out_dir  = 'C:\Users\hsollim\Documents\MATLAB\MATLAB\2p\extract_figures_truncated';

% Seconds (edit if needed)
skip_first_sec = 600;   % remove first 10 min
use_last_sec   = 600;   % then use last 10 min of what remains

run_id = sprintf('CS1014-1b · mirror2path · skip%ds_last%ds', skip_first_sec, use_last_sec);

params = struct( ...
    'z_thresh', 2, ...
    'n_boot', 1000, ...
    'time_bin', 1, ...
    'n_pc', 10, ...
    'n_corr_sample', 500, ...
    'frame_px', 2304, ...
    'max_show', 2000, ...
    'skip_first_sec', skip_first_sec, ...
    'use_last_sec', use_last_sec);

extract_analysis2_core(bl_path, inj_path, out_dir, run_id, params);

fprintf('\nTotal wall time (script): %.1f s\n', toc);
