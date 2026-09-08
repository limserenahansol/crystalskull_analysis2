function B = extract_bl_inj_neuron_bundle(bl_path, inj_path, params)
%EXTRACT_BL_INJ_NEURON_BUNDLE  Load BL/INJ traces, same crops as extract_analysis2_core, return per-neuron + pop data.
%   B = extract_bl_inj_neuron_bundle(bl_path, inj_path, params)
%
%   Used for cohort-level figures pooling many neurons across mice. Fields include
%   z_bl, z_inj, spatial weights for injection, population bin vectors, etc.

if nargin < 3 || isempty(params)
    params = struct();
end

z_thresh = bundle_getp(params, 'z_thresh', 2);
time_bin = bundle_getp(params, 'time_bin', 1);
frame_px = bundle_getp(params, 'frame_px', 2304);

Rb = extract_resolve_extract_files(bl_path);
Ri = extract_resolve_extract_files(inj_path);
if isempty(Rb.frr) || isempty(Rb.summary) || isempty(Ri.frr) || isempty(Ri.summary)
    error('extract_bl_inj_neuron_bundle:MissingFiles', 'Could not resolve mats under BL or INJ path.');
end

T_bl = double(h5read(Rb.frr, '/T'))';
S_inj_ir   = double(h5read(Ri.frr, '/S/ir'));
S_inj_jc   = double(h5read(Ri.frr, '/S/jc'));
S_inj_data = double(h5read(Ri.frr, '/S/data'));

ext_inj_file = Ri.ds_ext;
sum_bl_file  = Rb.summary;
sum_inj_file = Ri.summary;

T_inj = double(h5read(Ri.frr, '/T'))';
sumimg_inj = h5read(ext_inj_file, '/info/summary_image')';

fs_bl  = h5read(sum_bl_file,  '/framerate');
fs_inj = h5read(sum_inj_file, '/framerate');
fs_bl = double(fs_bl(1));
fs_inj = double(fs_inj(1));

[n_frames_bl,  n_neurons_bl]  = size(T_bl);
[n_frames_inj, n_neurons_inj] = size(T_inj);
dur_bl  = n_frames_bl  / fs_bl;
dur_inj = n_frames_inj / fs_inj;

[sk_bl, ul_bl, sk_inj, ul_inj, uf_bl, uf_inj] = extract_resolve_time_window_params(params);
match_inj = bundle_getp(params, 'match_inj_last_to_baseline_dur', false);

crop_bl  = (uf_bl > 0) || (sk_bl > 0) || (isfinite(ul_bl) && ul_bl > 0);
if crop_bl
    T_bl = extract_apply_time_window_T(T_bl, fs_bl, sk_bl, ul_bl, uf_bl);
end
[n_frames_bl, n_neurons_bl] = size(T_bl);
dur_bl = n_frames_bl / fs_bl;

if match_inj
    ul_inj_eff = dur_bl;
    uf_inj_eff = 0;
    crop_inj = (sk_inj > 0) || (isfinite(ul_inj_eff) && ul_inj_eff > 0);
    if crop_inj
        T_inj = extract_apply_time_window_T(T_inj, fs_inj, sk_inj, ul_inj_eff, uf_inj_eff);
    end
else
    crop_inj = (uf_inj > 0) || (sk_inj > 0) || (isfinite(ul_inj) && ul_inj > 0);
    if crop_inj
        T_inj = extract_apply_time_window_T(T_inj, fs_inj, sk_inj, ul_inj, uf_inj);
    end
end
[n_frames_inj, n_neurons_inj] = size(T_inj);
dur_inj = n_frames_inj / fs_inj;

nf = n_frames_bl;
win = max(round(nf * 0.15), 300);
half = floor(win / 2);
F0 = nan(size(T_bl));
for c = 1:half:nf
    s = max(1, c - half);
    e = min(nf, c + half);
    f0_chunk = prctile(T_bl(s:e, :), 8, 1);
    fe = min(c + half - 1, nf);
    F0(c:fe, :) = repmat(f0_chunk, fe - c + 1, 1);
end
for j = 1:size(F0, 2)
    col = F0(:, j);
    m = isnan(col);
    if any(m)
        iv = find(~m);
        col(m) = interp1(iv, col(iv), find(m), 'nearest', 'extrap');
        F0(:, j) = col;
    end
end
F0(F0 < 0.5) = 0.5;
dff_bl = (T_bl - F0) ./ F0;
clear F0 T_bl

nf = n_frames_inj;
win = max(round(nf * 0.15), 300);
half = floor(win / 2);
F0 = nan(size(T_inj));
for c = 1:half:nf
    s = max(1, c - half);
    e = min(nf, c + half);
    f0_chunk = prctile(T_inj(s:e, :), 8, 1);
    fe = min(c + half - 1, nf);
    F0(c:fe, :) = repmat(f0_chunk, fe - c + 1, 1);
end
for j = 1:size(F0, 2)
    col = F0(:, j);
    m = isnan(col);
    if any(m)
        iv = find(~m);
        col(m) = interp1(iv, col(iv), find(m), 'nearest', 'extrap');
        F0(:, j) = col;
    end
end
F0(F0 < 0.5) = 0.5;
dff_inj = (T_inj - F0) ./ F0;
clear F0 T_inj

z_bl  = (dff_bl  - mean(dff_bl, 1))  ./ max(std(dff_bl, 0, 1),  1e-8);
z_inj = (dff_inj - mean(dff_inj, 1)) ./ max(std(dff_inj, 0, 1), 1e-8);

mean_dff_bl   = mean(dff_bl, 1)';
mean_dff_inj  = mean(dff_inj, 1)';
std_dff_bl    = std(dff_bl, 0, 1)';
std_dff_inj   = std(dff_inj, 0, 1)';
skew_dff_bl   = skewness(dff_bl, 1, 1)';
skew_dff_inj  = skewness(dff_inj, 1, 1)';
clear dff_bl dff_inj

cross_bl  = diff(double(z_bl  > z_thresh), 1, 1);
cross_inj = diff(double(z_inj > z_thresh), 1, 1);
event_count_bl  = sum(cross_bl  == 1, 1)';
event_count_inj = sum(cross_inj == 1, 1)';
event_rate_bl   = event_count_bl  / dur_bl;
event_rate_inj  = event_count_inj / dur_inj;
clear cross_bl cross_inj

bin_frames_bl  = round(fs_bl  * time_bin);
bin_frames_inj = round(fs_inj * time_bin);
n_bins_bl  = floor(n_frames_bl  / bin_frames_bl);
n_bins_inj = floor(n_frames_inj / bin_frames_inj);

pop_rate_bl  = zeros(n_bins_bl, 1);
pop_frac_bl  = zeros(n_bins_bl, 1);
pop_rate_inj = zeros(n_bins_inj, 1);
pop_frac_inj = zeros(n_bins_inj, 1);

for b = 1:n_bins_bl
    idx = (b - 1) * bin_frames_bl + 1 : b * bin_frames_bl;
    chunk = z_bl(idx, :);
    pop_rate_bl(b) = mean(chunk(:));
    pop_frac_bl(b) = mean(any(chunk > z_thresh, 1));
end
for b = 1:n_bins_inj
    idx = (b - 1) * bin_frames_inj + 1 : b * bin_frames_inj;
    chunk = z_inj(idx, :);
    pop_rate_inj(b) = mean(chunk(:));
    pop_frac_inj(b) = mean(any(chunk > z_thresh, 1));
end

B = struct();
B.mean_dff_bl = mean_dff_bl;
B.mean_dff_inj = mean_dff_inj;
B.std_dff_bl = std_dff_bl;
B.std_dff_inj = std_dff_inj;
B.skew_dff_bl = skew_dff_bl;
B.skew_dff_inj = skew_dff_inj;
B.event_rate_bl = event_rate_bl;
B.event_rate_inj = event_rate_inj;
B.z_bl = z_bl;
B.z_inj = z_inj;
B.pop_rate_bl = pop_rate_bl;
B.pop_frac_bl = pop_frac_bl;
B.pop_rate_inj = pop_rate_inj;
B.pop_frac_inj = pop_frac_inj;
B.n_bins_bl = n_bins_bl;
B.n_bins_inj = n_bins_inj;
B.fs_bl = fs_bl;
B.fs_inj = fs_inj;
B.dur_bl = dur_bl;
B.dur_inj = dur_inj;
B.n_neurons_bl = n_neurons_bl;
B.n_neurons_inj = n_neurons_inj;
B.S_inj_ir = S_inj_ir;
B.S_inj_jc = S_inj_jc;
B.S_inj_data = S_inj_data;
B.frame_px = frame_px;
B.sumimg_inj = sumimg_inj;
B.z_thresh = z_thresh;
B.time_bin = time_bin;

end

function v = bundle_getp(s, name, default)
if isstruct(s) && isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
