function m = extract_session_metrics(session_mirror_path, params)
%EXTRACT_SESSION_METRICS  Population metrics for one recording (one mirror folder).
%
%   m = extract_session_metrics(session_mirror_path, params)
%
%   session_mirror_path — folder containing resolved *moco_frr.mat, etc. (e.g.
%     .../Baseline/mirror1path or .../Injection/mirror1path)
%
%   params — same optional fields as extract_analysis2_core (z_thresh, n_boot,
%     time_bin, n_pc, n_corr_sample, skip_first_sec, use_last_sec, frame_px).
%
%   Output m fields are scalars suitable for Δ = m_inj - m_baseline per mouse.

if nargin < 2 || isempty(params)
    params = struct();
end

z_thresh = getp(params, 'z_thresh', 2);
time_bin = getp(params, 'time_bin', 1);
n_pc     = getp(params, 'n_pc', 10);
n_corr_sample = getp(params, 'n_corr_sample', 500);
skip_first_sec = getp(params, 'skip_first_sec', 0);
use_last_sec   = getp(params, 'use_last_sec', inf);

R = extract_resolve_extract_files(session_mirror_path);
if isempty(R.frr) || isempty(R.ds_ext) || isempty(R.summary)
    error('extract_session_metrics:MissingFiles', ...
        'Could not resolve moco_frr / ds_ext / summary under:\n%s', session_mirror_path);
end

frr_file = R.frr;
sum_file = R.summary;

T = double(h5read(frr_file, '/T'))';

fs = h5read(sum_file, '/framerate');
fs = double(fs(:));
fs = fs(1);

[n_frames, n_neurons] = size(T);
dur = n_frames / fs;

if skip_first_sec > 0 || (isfinite(use_last_sec) && use_last_sec > 0)
    T = extract_apply_time_window_T(T, fs, skip_first_sec, use_last_sec);
    [n_frames, n_neurons] = size(T);
    dur = n_frames / fs;
end

% dF/F (8th pct rolling baseline)
nf = n_frames;
win = max(round(nf * 0.15), 300);
half = floor(win / 2);
F0 = nan(size(T));
for c = 1:half:nf
    s = max(1, c - half);
    e = min(nf, c + half);
    f0_chunk = prctile(T(s:e, :), 8, 1);
    fe = min(c + half - 1, nf);
    F0(c:fe, :) = repmat(f0_chunk, fe - c + 1, 1);
end
for jj = 1:size(F0, 2)
    col = F0(:, jj);
    mm = isnan(col);
    if any(mm)
        iv = find(~mm);
        col(mm) = interp1(iv, col(iv), find(mm), 'nearest', 'extrap');
        F0(:, jj) = col;
    end
end
F0(F0 < 0.5) = 0.5;
dff = (T - F0) ./ F0;
clear F0 T

z = (dff - mean(dff, 1)) ./ max(std(dff, 0, 1), 1e-8);

mean_dff_n = mean(dff, 1)';
std_dff_n  = std(dff, 0, 1)';
skew_dff_n = skewness(dff, 1, 1)';

cross_z = diff(double(z > z_thresh), 1, 1);
event_count = sum(cross_z == 1, 1)';
event_rate_n = event_count / dur;
clear cross_z

bin_frames = round(fs * time_bin);
n_bins = floor(n_frames / bin_frames);
pop_rate = zeros(n_bins, 1);
pop_frac = zeros(n_bins, 1);
for b = 1:n_bins
    idx = (b - 1) * bin_frames + 1 : b * bin_frames;
    chunk = z(idx, :);
    pop_rate(b) = mean(chunk(:));
    pop_frac(b) = mean(any(chunk > z_thresh, 1));
end

rng(42);
idx_c = randperm(n_neurons, min(n_corr_sample, n_neurons));
C = corrcoef(z(:, idx_c));
mask_triu = triu(true(size(C, 1)), 1);
pw = C(mask_triu);

Z_bin = zeros(n_bins, n_neurons);
for b = 1:n_bins
    Z_bin(b, :) = mean(z((b - 1) * bin_frames + 1 : b * bin_frames, :), 1);
end
[~, ~, ~, ~, expl] = pca(Z_bin, 'NumComponents', n_pc);
clear Z_bin z dff

m = struct();
m.session_path = session_mirror_path;
m.n_neurons = n_neurons;
m.n_frames = n_frames;
m.fs = fs;
m.dur_s = dur;
m.mean_mean_dff = mean(mean_dff_n);
m.median_mean_dff = median(mean_dff_n);
m.mean_std_dff = mean(std_dff_n);
m.median_std_dff = median(std_dff_n);
m.mean_skew_dff = mean(skew_dff_n);
m.median_skew_dff = median(skew_dff_n);
m.mean_event_rate = mean(event_rate_n);
m.median_event_rate = median(event_rate_n);
m.mean_pop_frac = mean(pop_frac);
m.mean_pop_rate_z = mean(pop_rate);
m.median_pw_corr = median(pw);
m.pc1_var_pct = expl(1);
m.pc1_to_pc3_var_pct = sum(expl(1:min(3, numel(expl))));

% Keep raw vectors for optional per-neuron Δ (same neuron count only — not used across days)
m.mean_dff_per_neuron = mean_dff_n;
m.event_rate_per_neuron = event_rate_n;

end

function v = getp(s, name, default)
if isstruct(s) && isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
