% EXTRACT Results: Baseline vs Injection — Population-Level Analysis (single mouse / mirror)
% Mesoscope / Crystal-skull whole dorsal cortex
%
% For multiple mice and mirrors, run extract_analysis2_batch.m instead.
%
% Requires: extract_analysis2_core.m on the MATLAB path.

clear; close all; clc;
tic;

%% --- Paths & parameters -------------------------------------------------
base_dir = '\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\CS1014-1b';
bl_path  = fullfile(base_dir, 'Baseline',  'mirror2path');
inj_path = fullfile(base_dir, 'Injection', 'mirror2path');
out_dir  = 'C:\Users\hsollim\Documents\MATLAB\MATLAB\2p\extract_figures';

% Label for figure titles and summary_stats (e.g. mouse ID, or mouse + mirror)
run_id = 'CS1014-1b · mirror2path';

params = struct( ...
    'z_thresh', 2, ...
    'n_boot', 1000, ...
    'time_bin', 1, ...
    'n_pc', 10, ...
    'n_corr_sample', 500, ...
    'frame_px', 2304, ...
    'max_show', 2000);

extract_analysis2_core(bl_path, inj_path, out_dir, run_id, params);

fprintf('\nTotal wall time (script): %.1f s\n', toc);
