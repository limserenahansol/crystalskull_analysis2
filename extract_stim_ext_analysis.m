function R = extract_stim_ext_analysis(ext_file, stim, out_dir, run_id, params)
%EXTRACT_STIM_EXT_ANALYSIS  Peri-poke analysis from an EXTRACT *_ext.mat file.
%
%   R = extract_stim_ext_analysis(ext_file, stim, out_dir, run_id, params)
%
% For the crystal-skull poke datasets where the cell-extraction output is a
% single EXTRACT *_moco_ds_pre_ext.mat (HDF5/v7.3) holding:
%     /temporal_weights  [frames x neurons]  (EXTRACT activity traces, dF/F-like)
%     /config/downsample_time_by             (temporal downsample factor)
%     /info/summary_image [px x px]          (FOV image)
% There is no separate moco_frr / summary .mat. Spatial footprints are
% stored as ndSparse; extract_load_ext_spatial reconstructs them from the
% HDF5 CSC under /#refs# so fig4 is a real poke-response map.
%
% INPUTS
%   ext_file : full path to the *_ext.mat file.
%   stim     : struct with movie-frame indices .peak/.pre/.post (vectors).
%   out_dir, run_id : output folder + label.
%   params   : optional struct. Extra fields beyond the core ones:
%                fps_movie   (33)    original movie frame rate (Hz)
%                recompute_dff (false) if true, apply 8th-pct rolling dF/F to
%                                    temporal_weights before z-scoring (usually
%                                    unnecessary: EXTRACT traces are dF/F-like)
%                event_field ('peak') 'peak' | 'pre' | 'post' -> t=0 anchor
%              plus core fields: peri_win, base_win, resp_win, z_thresh,
%              resp_test, alpha, max_show (see extract_stim_peristim_core.m).
%
% Requires extract_stim_peristim_core.m on the path.

if nargin < 5 || isempty(params); params = struct(); end
tic;

fps_movie     = getp(params, 'fps_movie', 33);
recompute_dff = getp(params, 'recompute_dff', false);
event_field   = getp(params, 'event_field', 'peak');
frame_px      = getp(params, 'frame_px', 2304);

if exist(ext_file, 'file') ~= 2
    error('extract_stim_ext_analysis:NoFile', 'EXTRACT file not found:\n%s', ext_file);
end

fprintf('=== STIM (ext) : %s ===\n  file: %s\n', run_id, ext_file);

% --- Trace matrix (frames x neurons) ------------------------------------
tw = h5read(ext_file, '/temporal_weights');
% Orient so rows = frames. Frames dimension must cover the max poke index.
maxidx = max(stim.peak(:));
if isfield(stim, 'pre') && ~isempty(stim.pre)
    maxidx = max(maxidx, max(stim.pre(:)));
end
if isfield(stim, 'post') && ~isempty(stim.post)
    maxidx = max(maxidx, max(stim.post(:)));
end
if size(tw, 1) < size(tw, 2) && size(tw, 2) >= maxidx && size(tw, 1) < maxidx
    tw = tw.';
elseif size(tw, 1) < maxidx && size(tw, 2) >= maxidx
    tw = tw.';
end
T = double(tw);
[n_frames, n_neurons] = size(T);

% --- Frame rate + onset mapping (movie frames -> trace frames) ----------
dsb = 1;
try
    dsb = double(h5read(ext_file, '/config/downsample_time_by'));
catch
end
if ~isfinite(dsb) || dsb <= 0; dsb = 1; end
fs = fps_movie / dsb;
fprintf('  %d neurons | %d frames | %.3f Hz (fps=%.1f / ds=%g)\n', ...
    n_neurons, n_frames, fs, fps_movie, dsb);

% --- Summary image (optional) -------------------------------------------
sumimg = [];
try
    sumimg = h5read(ext_file, '/info/summary_image')';
catch
end

% --- dF/F (optional) + z-score ------------------------------------------
if recompute_dff
    fprintf('  Recomputing dF/F (8th-pct rolling)...\n');
    win = max(round(n_frames*0.15), 300); half = floor(win/2);
    F0 = nan(size(T));
    for c = 1:half:n_frames
        s = max(1, c-half); e = min(n_frames, c+half);
        f0c = prctile(T(s:e, :), 8, 1);
        fe = min(c+half-1, n_frames);
        F0(c:fe, :) = repmat(f0c, fe-c+1, 1);
    end
    for j = 1:size(F0, 2)
        col = F0(:, j); m = isnan(col);
        if any(m); iv = find(~m); col(m) = interp1(iv, col(iv), find(m), 'nearest', 'extrap'); F0(:, j) = col; end
    end
    F0(F0 < 0.5) = 0.5;
    T = (T - F0) ./ F0;
    clear F0
end
z = (T - mean(T, 1)) ./ max(std(T, 0, 1), 1e-8);
clear T tw

% --- Resolve onset indices + pre/post offsets ---------------------------
switch lower(event_field)
    case 'peak'; anchor_raw = stim.peak;  ef_lbl = 'peak_idx';
    case 'pre';  anchor_raw = stim.pre;   ef_lbl = 'peak_pre_idx';
    case 'post'; anchor_raw = stim.post;  ef_lbl = 'peak_post_idx';
    otherwise; error('extract_stim_ext_analysis:BadEventField', 'event_field must be peak|pre|post.');
end
onset_movie = round(double(anchor_raw(:)));
onset_movie = onset_movie(isfinite(onset_movie) & onset_movie >= 1);
onset_trace = max(1, round(onset_movie / dsb));

% Markers: mean position of the pre- and post-minima relative to the anchor.
off = struct('ef_lbl', ef_lbl, ...
             'pre_s',  mean_offset_sec(stim.pre,  anchor_raw, fps_movie), ...
             'post_s', mean_offset_sec(stim.post, anchor_raw, fps_movie));

% --- Spatial footprints (ndSparse -> CSC under /#refs#) -----------------
Sn = [];
try
    [Sn, frameSize] = extract_load_ext_spatial(ext_file, n_neurons, frame_px);
    if ~isempty(frameSize)
        frame_px = frameSize(1);
    end
    fprintf('  Spatial footprints: %d x %d (normalized columns)\n', size(Sn, 1), size(Sn, 2));
catch ME
    warning('extract_stim_ext_analysis:NoSpatial', ...
        'Could not reconstruct spatial_weights (%s). Fig4 will be anatomy only.', ME.message);
end

R = extract_stim_peristim_core(z, fs, onset_trace, off, sumimg, Sn, frame_px, ...
    out_dir, run_id, params);
R.ext_file = ext_file;
save(fullfile(out_dir, 'stim_results.mat'), 'R', '-v7.3');
fprintf('Done (%s) in %.1f s -> %s\n', run_id, toc, out_dir);

end

function v = getp(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)); v = s.(name); else; v = default; end
end

function off = mean_offset_sec(target_list, anchor_list, fps)
%MEAN_OFFSET_SEC  Mean (target - anchor) in seconds, element-wise; NaN if unusable.
off = NaN;
if isempty(target_list) || isempty(anchor_list); return; end
a = double(target_list(:)); b = double(anchor_list(:));
if numel(a) ~= numel(b); return; end
off = mean((a - b) / fps);
end
