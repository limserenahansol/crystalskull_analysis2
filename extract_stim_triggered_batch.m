function extract_stim_triggered_batch()
%EXTRACT_STIM_TRIGGERED_BATCH  Peri-poke analysis for all movies in a stim file.
%
%   extract_stim_triggered_batch
%
% Loads Omer's estimated_stim_timing.mat, maps each movie (fnames{i}) to its
% EXTRACT output folder, runs extract_stim_triggered_analysis per movie, and
% writes a cohort summary (CSV + pooled population PSTH figure).
%
% estimated_stim_timing.mat is expected to contain:
%   fnames{i}        - movie file name (e.g. 'm201_d250813_..._mirror2.h5')
%   peak_idx{i}      - peak frame indices for movie i
%   peak_pre_idx{i}  - minimum frame indices just BEFORE each peak
%   peak_post_idx{i} - minimum frame indices just AFTER each peak
%
% Edit the CONFIG block, then run:  extract_stim_triggered_batch
%
% Requires on the path:
%   extract_stim_triggered_analysis.m, extract_resolve_extract_files.m

%% ========================= CONFIG ======================================
% Omer's stimulus-timing file:
stim_mat_file = '/mnt/crystal_skull/omer_hazon/estimated_stim_timing.mat';

% Parent folder(s) that contain the per-movie EXTRACT outputs. Each movie's
% moco_frr/ds_ext/summary .mat live in a subfolder somewhere under these roots.
ext_roots = { ...
    '/mnt/crystal_skull/omer_hazon/20250813', ...
    '/mnt/crystal_skull/omer_hazon/20250818'};

% Output root (one subfolder per movie is created here):
out_root = fullfile(pwd, 'stim_triggered_figures');

% Optional explicit map from movie file name -> EXTRACT folder. Leave empty
% to auto-discover by matching the movie basename against *moco_frr.mat files.
%   manual_map = containers.Map( ...
%     {'m201_d250813_s00cage4978415_m2c-fps33-mirror2.h5'}, ...
%     {'/mnt/crystal_skull/omer_hazon/20250813/m2c_mirror2'});
manual_map = containers.Map('KeyType', 'char', 'ValueType', 'char');

% Analysis parameters (see extract_stim_triggered_analysis.m header):
params = struct( ...
    'event_field', 'peak', ...   % 'peak' | 'pre' | 'post'  (t=0 anchor)
    'peri_win',  [-2 3], ...
    'base_win',  [-2 0], ...
    'resp_win',  [0 2], ...
    'z_thresh',  2, ...
    'resp_test', 'signrank', ...
    'alpha',     0.05, ...
    'frame_px',  2304, ...
    'fps_expected', 33, ...
    'use_single', true, ...
    'max_show',  2000);

% Limit to specific movie basenames (substring match); {} = all in fnames.
only_movies = {};
%% =======================================================================

close all; clc;

if exist(stim_mat_file, 'file') ~= 2
    error('extract_stim_triggered_batch:NoStimFile', ...
        'Stim timing file not found:\n%s', stim_mat_file);
end
if exist(out_root, 'dir') ~= 7; mkdir(out_root); end

fprintf('Loading stim timing: %s\n', stim_mat_file);
St = load_stim_timing(stim_mat_file);
n_mov = numel(St.fnames);
fprintf('  %d movie(s) in fnames.\n', n_mov);

results = {};
n_ok = 0; n_skip = 0; n_fail = 0;

for i = 1:n_mov
    fname = St.fnames{i};
    base  = strip_ext(fname);

    if ~isempty(only_movies) && ~any(cellfun(@(s) contains(base, s), only_movies))
        continue;
    end

    % Resolve EXTRACT folder for this movie
    if isKey(manual_map, fname)
        ext_folder = manual_map(fname);
    else
        ext_folder = find_ext_folder(ext_roots, base);
    end
    if isempty(ext_folder)
        fprintf('SKIP  %s  (no EXTRACT folder found)\n', base);
        n_skip = n_skip + 1;
        continue;
    end

    stim = struct('peak', St.peak_idx{i}, 'pre', St.peak_pre_idx{i}, ...
                  'post', St.peak_post_idx{i});
    out_dir = fullfile(out_root, base);
    run_id  = base;

    fprintf('\n-------- MOVIE %d/%d: %s --------\n', i, n_mov, base);
    try
        R = extract_stim_triggered_analysis(ext_folder, stim, out_dir, run_id, params);
        R.movie = base;
        results{end+1} = R; %#ok<AGROW>
        n_ok = n_ok + 1;
    catch ME
        warning('extract_stim_triggered_batch:RunFailed', 'FAILED %s: %s', base, ME.message);
        n_fail = n_fail + 1;
    end
    close all;
end

% ----- Cohort aggregation ------------------------------------------------
if ~isempty(results)
    write_cohort_summary(results, out_root, params);
end

fprintf('\n======== STIM BATCH DONE ========\n');
fprintf('  OK: %d | Skip: %d | Fail: %d\n', n_ok, n_skip, n_fail);
fprintf('  Outputs: %s\n', out_root);

end

% ======================= local helpers =================================
function St = load_stim_timing(matfile_path)
%LOAD_STIM_TIMING  Load fnames + peak index cells, normalizing types.
raw = load(matfile_path);
req = {'fnames', 'peak_idx', 'peak_pre_idx', 'peak_post_idx'};
for k = 1:numel(req)
    if ~isfield(raw, req{k})
        error('load_stim_timing:MissingVar', ...
            'Variable "%s" not found in %s', req{k}, matfile_path);
    end
end
St.fnames        = to_cellstr(raw.fnames);
St.peak_idx      = to_cellvec(raw.peak_idx);
St.peak_pre_idx  = to_cellvec(raw.peak_pre_idx);
St.peak_post_idx = to_cellvec(raw.peak_post_idx);
n = numel(St.fnames);
if numel(St.peak_idx) ~= n || numel(St.peak_pre_idx) ~= n || numel(St.peak_post_idx) ~= n
    error('load_stim_timing:LengthMismatch', ...
        'fnames (%d) and peak_* cell arrays must have equal length.', n);
end
end

function c = to_cellstr(x)
if iscell(x)
    c = cellfun(@char, x, 'UniformOutput', false);
elseif isstring(x)
    c = cellstr(x);
elseif ischar(x)
    c = cellstr(x);      % char matrix -> one row per name
else
    error('to_cellstr:BadType', 'fnames must be cell/string/char.');
end
c = c(:)';
end

function c = to_cellvec(x)
if iscell(x)
    c = cellfun(@(v) double(v(:)), x, 'UniformOutput', false);
    c = c(:)';
else
    c = {double(x(:))};   % single numeric vector -> one movie
end
end

function b = strip_ext(fname)
[~, b, ~] = fileparts(char(fname));
end

function folder = find_ext_folder(ext_roots, base)
%FIND_EXT_FOLDER  Locate the subfolder whose moco_frr .mat matches "base".
folder = '';
% Token from movie basename used for fuzzy matching (drop the h5 suffix bits).
key = regexprep(base, '\.h5$', '');
best_len = 0;
for r = 1:numel(ext_roots)
    root = ext_roots{r};
    if exist(root, 'dir') ~= 7; continue; end

    % Candidate folders: the root itself + all subfolders (2 levels deep).
    cands = collect_dirs(root, 2);
    for c = 1:numel(cands)
        S = extract_resolve_extract_files(cands{c});
        if isempty(S.frr); continue; end
        [~, frrname] = fileparts(S.frr);
        % Match if the movie key is a prefix of the frr file (or vice versa),
        % or the folder name contains the key.
        [~, foldername] = fileparts(cands{c});
        cand_keys = {frrname, foldername};
        for kk = 1:numel(cand_keys)
            ck = cand_keys{kk};
            L = longest_common_token(key, ck);
            if (startsWith(frrname, key) || contains(frrname, key) || ...
                    contains(foldername, key) || contains(key, ck)) && L > best_len
                folder = cands{c};
                best_len = L;
            end
        end
    end
end
end

function L = longest_common_token(a, b)
%LONGEST_COMMON_TOKEN  Heuristic overlap score between two names.
ta = split_tokens(a); tb = split_tokens(b);
common = intersect(ta, tb);
L = sum(cellfun(@numel, common));
end

function t = split_tokens(s)
t = regexp(s, '[_\-\.]', 'split');
t = t(~cellfun(@isempty, t));
end

function dirs = collect_dirs(root, depth)
dirs = {root};
if depth <= 0; return; end
d = dir(root);
for i = 1:numel(d)
    if ~d(i).isdir || strcmp(d(i).name, '.') || strcmp(d(i).name, '..'); continue; end
    sub = fullfile(root, d(i).name);
    dirs = [dirs, collect_dirs(sub, depth-1)]; %#ok<AGROW>
end
end

function write_cohort_summary(results, out_root, params)
%WRITE_COHORT_SUMMARY  CSV of per-movie scalars + pooled PSTH figure.
n = numel(results);

% CSV
csv_path = fullfile(out_root, 'cohort_stim_summary.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, ['movie,event_field,fs,n_neurons,n_pokes,n_trials,' ...
    'n_exc,frac_exc,n_inh,frac_inh,mean_amp,mean_pop_resp,p_pop,iti_median\n']);
for i = 1:n
    R = results{i};
    fprintf(fid, '%s,%s,%.4f,%d,%d,%d,%d,%.4f,%d,%.4f,%.5f,%.5f,%.3e,%.3f\n', ...
        R.movie, R.event_field, R.fs, R.n_neurons, R.n_pokes_total, R.n_trials, ...
        R.n_exc, R.frac_exc, R.n_inh, R.frac_inh, R.mean_amp, R.mean_pop_resp, ...
        R.p_pop, R.iti_median);
end
fclose(fid);
fprintf('Cohort CSV -> %s\n', csv_path);

% Pooled population PSTH on a common time grid
tq = params.peri_win(1) : 0.05 : params.peri_win(2);
psth = nan(n, numel(tq));
frac = nan(n, numel(tq));
for i = 1:n
    R = results{i};
    psth(i, :) = interp1(R.t_peri, R.pop_psth, tq, 'linear', NaN);
    frac(i, :) = interp1(R.t_peri, R.frac_psth, tq, 'linear', NaN);
end

figure('Position', [30 30 1500 560], 'Color', 'w');
sgtitle(sprintf('Cohort peri-poke response (%d movies, t=0 at %s)', ...
    n, results{1}.event_field), 'FontSize', 13, 'FontWeight', 'bold');

subplot(1, 2, 1);
plot(tq, psth', 'Color', [0.7 0.7 0.7 0.6], 'LineWidth', 0.8); hold on;
mu = mean(psth, 1, 'omitnan');
plot(tq, mu, 'r', 'LineWidth', 2.5);
xline(0, 'k-', 'LineWidth', 1.5);
xlabel('Time from poke (s)'); ylabel('Population mean z');
title('Per-movie (gray) + cohort mean (red)'); box off;

subplot(1, 2, 2);
plot(tq, frac', 'Color', [0.7 0.7 0.7 0.6], 'LineWidth', 0.8); hold on;
mu = mean(frac, 1, 'omitnan');
plot(tq, mu, 'k', 'LineWidth', 2.5);
xline(0, 'k-', 'LineWidth', 1.5);
xlabel('Time from poke (s)'); ylabel(sprintf('Fraction active (z>%g)', params.z_thresh));
title('Population recruitment'); box off;

exportgraphics(gcf, fullfile(out_root, 'cohort_stim_population.png'), 'Resolution', 200);
end
