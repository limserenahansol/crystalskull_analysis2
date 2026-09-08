function extract_stim_ext_batch(stim_mat_file, ext_root, out_root, params)
%EXTRACT_STIM_EXT_BATCH  Peri-poke analysis for all movies in a stim file (_ext layout).
%
%   extract_stim_ext_batch                              % uses CONFIG defaults below
%   extract_stim_ext_batch(stim_mat_file, ext_root, out_root, params)
%
% Loads Omer's estimated_stim_timing_*.mat, matches each movie (fnames{i}) to
% its EXTRACT *_ext.mat in ext_root, runs extract_stim_ext_analysis per movie,
% and writes a cohort summary (CSV + pooled population PSTH figure).
%
% estimated_stim_timing_*.mat must contain fnames + peak_idx + peak_pre_idx +
% peak_post_idx (one cell/entry per movie). fnames may end in "_moco.h5"; the
% matching extraction file is "<stem>*_ext.mat" in ext_root.
%
% Requires: extract_stim_ext_analysis.m, extract_stim_peristim_core.m.

%% ========================= CONFIG (used when called with no args) ======
if nargin < 1 || isempty(stim_mat_file)
    stim_mat_file = '\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\20250813\estimated_stim_timing_20250813.mat';
end
if nargin < 2 || isempty(ext_root)
    ext_root = fileparts(stim_mat_file);   % default: same folder as the stim file
end
if nargin < 3 || isempty(out_root)
    out_root = fullfile(fileparts(mfilename('fullpath')), 'stim_triggered_figures', ...
        ['ext_' datestr(now, 'yyyymmdd_HHMMSS')]);
end
if nargin < 4 || isempty(params)
    params = struct( ...
        'event_field', 'peak', ...   % 'peak' | 'pre' | 'post'  (t=0 anchor)
        'fps_movie',   33, ...
        'recompute_dff', false, ...  % EXTRACT temporal_weights already dF/F-like
        'peri_win',  [-2 3], ...
        'base_win',  [-2 0], ...
        'resp_win',  [0 2], ...
        'z_thresh',  2, ...
        'resp_test', 'signrank', ...
        'alpha',     0.05, ...
        'frame_px',  2304, ...
        'max_show',  2000);
end
only_movies = {};   % {} = all; else cell of basename substrings
%% =======================================================================

if exist(stim_mat_file, 'file') ~= 2
    error('extract_stim_ext_batch:NoStimFile', 'Stim file not found:\n%s', stim_mat_file);
end
if exist(out_root, 'dir') ~= 7; mkdir(out_root); end

fprintf('Stim file : %s\nExt root  : %s\nOut root  : %s\n', stim_mat_file, ext_root, out_root);
St = load_stim_timing(stim_mat_file);
n_mov = numel(St.fnames);
fprintf('  %d movie(s) in fnames.  event_field = %s\n', n_mov, params.event_field);

results = {}; n_ok = 0; n_skip = 0; n_fail = 0;
for i = 1:n_mov
    base = strip_ext(St.fnames{i});
    if ~isempty(only_movies) && ~any(cellfun(@(s) contains(base, s), only_movies)); continue; end

    ext_file = find_ext_file(ext_root, base);
    if isempty(ext_file)
        fprintf('SKIP  %s  (no *_ext.mat found)\n', base); n_skip = n_skip + 1; continue;
    end

    stim = struct('peak', St.peak_idx{i}, 'pre', St.peak_pre_idx{i}, 'post', St.peak_post_idx{i});
    out_dir = fullfile(out_root, base);

    fprintf('\n-------- MOVIE %d/%d: %s --------\n', i, n_mov, base);
    try
        R = extract_stim_ext_analysis(ext_file, stim, out_dir, base, params);
        R.movie = base; results{end+1} = R; %#ok<AGROW>
        n_ok = n_ok + 1;
    catch ME
        warning('extract_stim_ext_batch:RunFailed', 'FAILED %s: %s', base, ME.message);
        n_fail = n_fail + 1;
    end
    close all;
end

if ~isempty(results); write_cohort_summary(results, out_root, params); end

fprintf('\n======== STIM EXT BATCH DONE ========\n');
fprintf('  OK: %d | Skip: %d | Fail: %d\n  Outputs: %s\n', n_ok, n_skip, n_fail, out_root);

end

% ======================= local helpers =================================
function St = load_stim_timing(matfile_path)
raw = load(matfile_path);
if ~isfield(raw, 'fnames') || ~isfield(raw, 'peak_idx')
    error('load_stim_timing:MissingVar', 'Need fnames and peak_idx in %s', matfile_path);
end
St.fnames   = to_cellstr(raw.fnames);
St.peak_idx = to_cellvec(raw.peak_idx);
n = numel(St.fnames);
if numel(St.peak_idx) ~= n
    error('load_stim_timing:LengthMismatch', 'fnames and peak_idx must be equal length.');
end
% 20250818 has dip-based peak_idx only; 20250813 also has pre/post minima.
if isfield(raw, 'peak_pre_idx')
    St.peak_pre_idx = to_cellvec(raw.peak_pre_idx);
else
    St.peak_pre_idx = repmat({[]}, 1, n);
end
if isfield(raw, 'peak_post_idx')
    St.peak_post_idx = to_cellvec(raw.peak_post_idx);
else
    St.peak_post_idx = repmat({[]}, 1, n);
end
if numel(St.peak_pre_idx) ~= n
    St.peak_pre_idx = repmat({[]}, 1, n);
end
if numel(St.peak_post_idx) ~= n
    St.peak_post_idx = repmat({[]}, 1, n);
end
end

function c = to_cellstr(x)
if iscell(x); c = cellfun(@char, x, 'UniformOutput', false);
elseif isstring(x); c = cellstr(x);
elseif ischar(x); c = cellstr(x);
else; error('to_cellstr:BadType', 'fnames must be cell/string/char.'); end
c = c(:)';
end

function c = to_cellvec(x)
if iscell(x); c = cellfun(@(v) double(v(:)), x, 'UniformOutput', false); c = c(:)';
else; c = {double(x(:))}; end
end

function b = strip_ext(fname)
[~, b, ~] = fileparts(char(fname));
end

function ext_file = find_ext_file(ext_root, base)
%FIND_EXT_FILE  Locate "<base>*_ext.mat" in ext_root (exact-prefix preferred).
ext_file = '';
d = dir(fullfile(ext_root, '*_ext.mat'));
if isempty(d); return; end
names = {d.name};
% Prefer files whose name starts with the movie base stem.
hit = find(startsWith(names, base), 1);
if isempty(hit)
    % fall back: base contains a common token with the file stem
    hit = find(cellfun(@(n) contains(n, base) || contains(base, regexprep(n, '_ds_pre_ext\.mat$|_ext\.mat$', '')), names), 1);
end
if ~isempty(hit); ext_file = fullfile(ext_root, names{hit}); end
end

function write_cohort_summary(results, out_root, params)
n = numel(results);
csv_path = fullfile(out_root, 'cohort_stim_summary.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, ['movie,event_field,fs,n_neurons,n_pokes,n_trials,' ...
    'n_exc,frac_exc,n_inh,frac_inh,mean_amp,mean_pop_resp,p_pop,iti_median\n']);
for i = 1:n
    R = results{i};
    fprintf(fid, '%s,%s,%.4f,%d,%d,%d,%d,%.4f,%d,%.4f,%.5f,%.5f,%.3e,%.3f\n', ...
        R.movie, R.event_field, R.fs, R.n_neurons, R.n_pokes_total, R.n_trials, ...
        R.n_exc, R.frac_exc, R.n_inh, R.frac_inh, R.mean_amp, R.mean_pop_resp, R.p_pop, R.iti_median);
end
fclose(fid);
fprintf('Cohort CSV -> %s\n', csv_path);

tq = params.peri_win(1) : 0.05 : params.peri_win(2);
psth = nan(n, numel(tq)); frac = nan(n, numel(tq));
for i = 1:n
    R = results{i};
    psth(i, :) = interp1(R.t_peri, R.pop_psth, tq, 'linear', NaN);
    frac(i, :) = interp1(R.t_peri, R.frac_psth, tq, 'linear', NaN);
end

figure('Position', [30 30 1500 560], 'Color', 'w');
sgtitle(sprintf('Cohort peri-poke response (%d movies, t=0 at %s)', n, results{1}.event_field), ...
    'FontSize', 13, 'FontWeight', 'bold');
subplot(1, 2, 1);
plot(tq, psth', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); hold on;
plot(tq, mean(psth, 1, 'omitnan'), 'r', 'LineWidth', 2.5);
xline(0, 'k-', 'LineWidth', 1.5);
xlabel('Time from poke (s)'); ylabel('Population mean z');
title('Per-movie (gray) + cohort mean (red)'); box off;
subplot(1, 2, 2);
plot(tq, frac', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); hold on;
plot(tq, mean(frac, 1, 'omitnan'), 'k', 'LineWidth', 2.5);
xline(0, 'k-', 'LineWidth', 1.5);
xlabel('Time from poke (s)'); ylabel(sprintf('Fraction active (z>%g)', params.z_thresh));
title('Population recruitment'); box off;
exportgraphics(gcf, fullfile(out_root, 'cohort_stim_population.png'), 'Resolution', 200);
extract_stim_lab_meeting_figures(results, out_root, params);
end
