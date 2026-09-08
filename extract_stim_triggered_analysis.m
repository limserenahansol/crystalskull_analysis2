function R = extract_stim_triggered_analysis(ext_folder, stim, out_dir, run_id, params)
%EXTRACT_STIM_TRIGGERED_ANALYSIS  Peri-poke analysis from a moco_frr session folder.
%
%   R = extract_stim_triggered_analysis(ext_folder, stim, out_dir, run_id, params)
%
% Loader for the EXTRACT "moco_frr" layout (full-rate traces /T in
% *moco_frr.mat + sparse footprints /S + framerate in *summary.mat). Computes
% dF/F (8th-pct rolling baseline) and z-scores, then delegates the
% peri-stimulus math + figures to extract_stim_peristim_core.m (which also
% draws the spatial response map from the footprints).
%
% For the crystal-skull poke datasets whose extraction output is instead a
% single *_ext.mat with /temporal_weights, use extract_stim_ext_analysis.m.
%
% INPUTS
%   ext_folder : folder with moco_frr / ds_ext / summary .mat (resolved via
%                extract_resolve_extract_files.m).
%   stim       : struct with movie-frame indices .peak / .pre / .post.
%   out_dir, run_id : output folder + label.
%   params     : optional struct. Fields:
%                  event_field ('peak') 'peak'|'pre'|'post' -> t=0 anchor
%                  fps_expected (33)    sanity-check warning only
%                  n_frames_movie ([])  if set and != trace frames, rescale
%                                       stim indices by trace/movie ratio
%                  use_single (true)    store traces in single precision
%                  frame_px (2304)      FOV size for spatial map
%                plus core fields: peri_win, base_win, resp_win, z_thresh,
%                resp_test, alpha, max_show (see extract_stim_peristim_core.m).
%
% Requires: extract_resolve_extract_files.m, extract_stim_peristim_core.m.

if nargin < 5 || isempty(params); params = struct(); end
tic;

event_field  = getp(params, 'event_field', 'peak');
fps_expected = getp(params, 'fps_expected', 33);
use_single   = getp(params, 'use_single', true);
frame_px     = getp(params, 'frame_px', 2304);

Rf = extract_resolve_extract_files(ext_folder);
if isempty(Rf.frr) || isempty(Rf.summary)
    error('extract_stim_triggered_analysis:MissingFiles', ...
        'Could not find moco_frr / summary .mat under:\n%s', ext_folder);
end
fprintf('=== STIM (frr): %s ===\n  frr: %s\n', run_id, Rf.frr);

T = h5read(Rf.frr, '/T')';                 % frames x neurons
S_ir = double(h5read(Rf.frr, '/S/ir'));
S_jc = double(h5read(Rf.frr, '/S/jc'));
S_data = double(h5read(Rf.frr, '/S/data'));
fs = double(h5read(Rf.summary, '/framerate'));
if use_single; T = single(T); end
[n_frames, n_neurons] = size(T);
if abs(fs - fps_expected) > 1
    warning('extract_stim_triggered_analysis:FpsMismatch', ...
        'Trace framerate %.3f Hz differs from expected %.1f Hz.', fs, fps_expected);
end

sumimg = [];
if ~isempty(Rf.ds_ext)
    try
        sumimg = h5read(Rf.ds_ext, '/info/summary_image')';
    catch
        sumimg = [];
    end
end

% --- dF/F (8th-pct rolling baseline) + z-score --------------------------
fprintf('  Computing dF/F + z-score...\n');
win = max(round(n_frames*0.15), 300); half = floor(win/2);
F0 = nan(size(T), 'like', T);
for c = 1:half:n_frames
    s = max(1, c-half); e = min(n_frames, c+half);
    f0c = prctile(T(s:e, :), 8, 1);
    fe = min(c+half-1, n_frames);
    F0(c:fe, :) = repmat(cast(f0c, 'like', T), fe-c+1, 1);
end
for j = 1:size(F0, 2)
    col = F0(:, j); m = isnan(col);
    if any(m); iv = find(~m); col(m) = interp1(iv, col(iv), find(m), 'nearest', 'extrap'); F0(:, j) = col; end
end
F0(F0 < 0.5) = 0.5;
dff = (T - F0) ./ F0;
clear F0 T
z = (dff - mean(dff, 1)) ./ max(std(dff, 0, 1), 1e-8);
clear dff

% --- Build normalized sparse footprints for the spatial map -------------
n_pix = frame_px * frame_px;
col_idx = zeros(S_jc(end), 1);
for c = 1:n_neurons; col_idx(S_jc(c)+1:S_jc(c+1)) = c; end
Ss = sparse(S_ir+1, col_idx, S_data, n_pix, n_neurons);
cs = full(sum(Ss, 1)); cs(cs == 0) = 1; Sn = Ss ./ cs;
clear Ss col_idx S_ir S_jc S_data

% --- Resolve onset indices (+ optional movie->trace rescale) ------------
switch lower(event_field)
    case 'peak'; anchor_raw = stim.peak;  ef_lbl = 'peak_idx';
    case 'pre';  anchor_raw = stim.pre;   ef_lbl = 'peak_pre_idx';
    case 'post'; anchor_raw = stim.post;  ef_lbl = 'peak_post_idx';
    otherwise; error('extract_stim_triggered_analysis:BadEventField', 'event_field must be peak|pre|post.');
end
onset = round(double(anchor_raw(:)));
onset = onset(isfinite(onset) & onset >= 1);
n_frames_movie = getp(params, 'n_frames_movie', []);
if ~isempty(n_frames_movie) && n_frames_movie > 0 && n_frames_movie ~= n_frames
    sc = n_frames / n_frames_movie;
    onset = max(1, round(onset * sc));
    fprintf('  Scaled stim indices by %.4f (movie %d -> trace %d).\n', sc, n_frames_movie, n_frames);
end

off = struct('ef_lbl', ef_lbl, ...
             'pre_s',  mean_offset_sec(stim.pre,  anchor_raw, fps_expected), ...
             'post_s', mean_offset_sec(stim.post, anchor_raw, fps_expected));

R = extract_stim_peristim_core(z, fs, onset, off, sumimg, Sn, frame_px, ...
    out_dir, run_id, params);
R.ext_folder = ext_folder;
fprintf('Done (%s) in %.1f s -> %s\n', run_id, toc, out_dir);

end

function v = getp(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)); v = s.(name); else; v = default; end
end

function off = mean_offset_sec(target_list, anchor_list, fps)
off = NaN;
if isempty(target_list) || isempty(anchor_list); return; end
a = double(target_list(:)); b = double(anchor_list(:));
if numel(a) ~= numel(b); return; end
off = mean((a - b) / fps);
end
