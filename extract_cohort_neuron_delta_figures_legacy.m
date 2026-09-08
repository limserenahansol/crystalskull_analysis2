function extract_cohort_neuron_delta_figures_legacy()
%EXTRACT_COHORT_NEURON_DELTA_FIGURES_LEGACY  Reproduction of older cohort-neuron Δ export (pre output_granularity).
%
%   Use this to match archived outputs (e.g. cohort_neuron_delta_stats.txt with n_mice=, no Time/Scope lines,
%   population title "Population Δ (inj - baseline session mean)", SEM "across mice").
%
%   Behavior (matches pre–three-mode layout):
%     compare_mirrors_separately = true  → one output folder per mirror name; both mice pooled in each folder.
%     compare_mirrors_separately = false → single out_dir; all mouse×mirror bundles pooled.
%
%   Same math as extract_cohort_neuron_delta_figures.m (bundle load, Δ definitions, figures).
%   Requires: extract_bl_inj_neuron_bundle.m, extract_cohort_normalize_mirror_list.m, extract_resolve_* , etc.

tic;

%% ========================= CONFIG (edit like the main script) ==================
cohort_a_name = 'melittin_only';
cohort_a_root = '\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon';
cohort_a_mice = {'CS1014-1b', 'CS1014-2w'};

cohort_b_name = 'melittin_morphine';
cohort_b_root = '\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\morphine_20260408';
cohort_b_mice = {'1014m01b', '1014m02w'};

baseline_folder = 'Baseline';
injection_folder = 'Injection';
mirror_subfolders = {'mirror1path', 'mirror2path'};

% true  = out_dir/<mirror>/... with both mice pooled (typical old archive).
% false = all mirrors × mice in one out_dir.
compare_mirrors_separately = true;

time_window_mode = 'matched';  % set 'whole' or 'matched' to match the run you compare

out_dir_parent = 'C:\Users\hsollim\Documents\MATLAB\MATLAB\2p\cohort_neuron_delta_melittin_vs_morphine_legacy_check';
out_dir_override = '';
if isempty(out_dir_override)
    out_dir = fullfile(out_dir_parent, time_window_mode);
else
    out_dir = out_dir_override;
end

params = cohort_neuron_figure_params_legacy(time_window_mode);

col_a = [0.13 0.59 0.95];
col_b = [1.00 0.34 0.13];
max_show = 2000;
%% ==============================================================================

close all;
clc;

if strcmpi(time_window_mode, 'matched')
    tw_label = 'Time: INJ last = baseline duration (matched)';
else
    tw_label = 'Time: whole BL + whole INJ';
end
fprintf('extract_cohort_neuron_delta_figures_legacy | mode=%s | %s\n', time_window_mode, tw_label);
fprintf('Base output: %s\n', out_dir);

mirror_list = extract_cohort_normalize_mirror_list(mirror_subfolders);
n_mice_a = numel(cohort_a_mice);
n_mice_b = numel(cohort_b_mice);

if compare_mirrors_separately
    for km = 1:numel(mirror_list)
        mname = mirror_list{km};
        fprintf('\n--- Focal plane: %s ---\n', mname);
        fprintf('Loading cohort A...\n');
        [Ba, ~] = collect_bundles_legacy(cohort_a_root, cohort_a_mice, baseline_folder, injection_folder, mname, params);
        fprintf('Loading cohort B...\n');
        [Bb, ~] = collect_bundles_legacy(cohort_b_root, cohort_b_mice, baseline_folder, injection_folder, mname, params);
        if isempty(Ba) || isempty(Bb)
            warning('extract_cohort_neuron_delta_figures_legacy:SkipMirror', 'Skip mirror %s (empty cohort A or B).', mname);
            continue;
        end
        out_m = fullfile(out_dir, mname);
        if exist(out_m, 'dir') ~= 7
            mkdir(out_m);
        end
        cohort_neuron_delta_export_legacy(out_m, Ba, Bb, cohort_a_name, cohort_b_name, ...
            params, col_a, col_b, max_show, tw_label, n_mice_a, n_mice_b);
    end
else
    if exist(out_dir, 'dir') ~= 7
        mkdir(out_dir);
    end
    Ba = {};
    Bb = {};
    for km = 1:numel(mirror_list)
        mname = mirror_list{km};
        fprintf('Loading cohort A | %s...\n', mname);
        [b1, ~] = collect_bundles_legacy(cohort_a_root, cohort_a_mice, baseline_folder, injection_folder, mname, params);
        Ba = [Ba, b1]; %#ok<AGROW>
        fprintf('Loading cohort B | %s...\n', mname);
        [b2, ~] = collect_bundles_legacy(cohort_b_root, cohort_b_mice, baseline_folder, injection_folder, mname, params);
        Bb = [Bb, b2]; %#ok<AGROW>
    end
    if isempty(Ba) || isempty(Bb)
        error('extract_cohort_neuron_delta_figures_legacy:NoData', 'Need at least one session bundle per cohort.');
    end
    cohort_neuron_delta_export_legacy(out_dir, Ba, Bb, cohort_a_name, cohort_b_name, ...
        params, col_a, col_b, max_show, tw_label, n_mice_a, n_mice_b);
end

fprintf('\nDone (legacy). Elapsed %.1f s\n', toc);
end

function cohort_neuron_delta_export_legacy(out_dir, Ba, Bb, cohort_a_name, cohort_b_name, ...
    params, col_a, col_b, max_show, tw_label, n_mice_a, n_mice_b)

[d_mdf_a, d_er_a, d_std_a, d_sk_a, d_mdf_b, d_er_b, d_std_b, d_sk_b] = deal([]);
for k = 1:numel(Ba)
    B = Ba{k};
    d_mdf_a = [d_mdf_a; B.mean_dff_inj - mean(B.mean_dff_bl)]; %#ok<AGROW>
    d_er_a  = [d_er_a;  B.event_rate_inj - mean(B.event_rate_bl)];
    d_std_a = [d_std_a; B.std_dff_inj - mean(B.std_dff_bl)];
    d_sk_a  = [d_sk_a;  B.skew_dff_inj - mean(B.skew_dff_bl)];
end
for k = 1:numel(Bb)
    B = Bb{k};
    d_mdf_b = [d_mdf_b; B.mean_dff_inj - mean(B.mean_dff_bl)]; %#ok<AGROW>
    d_er_b  = [d_er_b;  B.event_rate_inj - mean(B.event_rate_bl)];
    d_std_b = [d_std_b; B.std_dff_inj - mean(B.std_dff_bl)];
    d_sk_b  = [d_sk_b;  B.skew_dff_inj - mean(B.skew_dff_bl)];
end

time_bin = params.time_bin;
[all_dpf_a, all_dpf_b, all_dpr_a, all_dpr_b, t_vec, mu_pf_a, sem_pf_a, mu_pf_b, sem_pf_b, mu_pr_a, sem_pr_a, mu_pr_b, sem_pr_b] = ...
    population_delta_series_legacy(Ba, Bb, time_bin);

fid = fopen(fullfile(out_dir, 'cohort_neuron_delta_stats.txt'), 'w');
fprintf(fid, 'Cohort neuron-level Δ: injection neuron metric minus MEAN of baseline session (same-cell pairing not required).\n');
fprintf(fid, 'Cohort A: %s | n_mice=%d | n_neurons_pooled=%d\n', cohort_a_name, n_mice_a, numel(d_mdf_a));
fprintf(fid, 'Cohort B: %s | n_mice=%d | n_neurons_pooled=%d\n', cohort_b_name, n_mice_b, numel(d_mdf_b));
fprintf(fid, 'CAVEAT: neurons cluster within mice; KS/ranksum treat each neuron as i.i.d. — interpret as exploratory or use mixed models for inference.\n\n');

log_cmp_legacy(fid, 'd_mean_dF/F', d_mdf_a, d_mdf_b);
log_cmp_legacy(fid, 'd_event_rate', d_er_a, d_er_b);
log_cmp_legacy(fid, 'd_std_dF/F', d_std_a, d_std_b);
log_cmp_legacy(fid, 'd_skewness', d_sk_a, d_sk_b);
log_cmp_legacy(fid, 'd_pop_frac (all 1s bins pooled)', all_dpf_a, all_dpf_b);
fclose(fid);

figure('Color', 'w', 'Position', [30 30 1500 900]);
sgtitle({sprintf(['Population Δ (inj - baseline session mean) | %s (n=%d neu) vs %s (n=%d neu)'], ...
    cohort_a_name, numel(d_mdf_a), cohort_b_name, numel(d_mdf_b)), tw_label}, 'FontSize', 11, 'FontWeight', 'bold');

subplot(2, 2, 1);
if ~isempty(t_vec)
    hold on;
    fill([t_vec; flipud(t_vec)], [mu_pr_a(:) - sem_pr_a(:); flipud(mu_pr_a(:) + sem_pr_a(:))], col_a, 'FaceAlpha', 0.22, 'EdgeColor', 'none');
    fill([t_vec; flipud(t_vec)], [mu_pr_b(:) - sem_pr_b(:); flipud(mu_pr_b(:) + sem_pr_b(:))], col_b, 'FaceAlpha', 0.22, 'EdgeColor', 'none');
    h1 = plot(t_vec, mu_pr_a, 'Color', col_a, 'LineWidth', 1.8);
    h2 = plot(t_vec, mu_pr_b, 'Color', col_b, 'LineWidth', 1.8);
    legend([h1, h2], {cohort_a_name, cohort_b_name}, 'Location', 'best');
end
xlabel('Time (s)');
ylabel('Δ pop. mean z');
title('Mean z-score: time courses (mean ± SEM across mice)');
box off;

subplot(2, 2, 2);
if ~isempty(t_vec)
    hold on;
    fill([t_vec; flipud(t_vec)], [mu_pf_a(:) - sem_pf_a(:); flipud(mu_pf_a(:) + sem_pf_a(:))], col_a, 'FaceAlpha', 0.22, 'EdgeColor', 'none');
    fill([t_vec; flipud(t_vec)], [mu_pf_b(:) - sem_pf_b(:); flipud(mu_pf_b(:) + sem_pf_b(:))], col_b, 'FaceAlpha', 0.22, 'EdgeColor', 'none');
    h1 = plot(t_vec, mu_pf_a, 'Color', col_a, 'LineWidth', 1.8);
    h2 = plot(t_vec, mu_pf_b, 'Color', col_b, 'LineWidth', 1.8);
    legend([h1, h2], {cohort_a_name, cohort_b_name}, 'Location', 'best');
end
xlabel('Time (s)');
ylabel(sprintf('Δ fraction active (z>%g)', params.z_thresh));
title('Recruitment: Δ time courses (mean ± SEM across mice)');
box off;

subplot(2, 2, 3);
edges = linspace(prctile([all_dpf_a; all_dpf_b], 1), prctile([all_dpf_a; all_dpf_b], 99), 50);
histogram(all_dpf_a, edges, 'Normalization', 'pdf', 'FaceColor', col_a, 'FaceAlpha', 0.45, 'EdgeColor', 'none');
hold on;
histogram(all_dpf_b, edges, 'Normalization', 'pdf', 'FaceColor', col_b, 'FaceAlpha', 0.45, 'EdgeColor', 'none');
xlabel('Δ synchrony (per 1 s bin)');
ylabel('Density');
[p_ks, ~] = kstest2(all_dpf_a, all_dpf_b);
[p_rs, ~] = ranksum(all_dpf_a, all_dpf_b);
title(sprintf('Pooled Δ synchrony bins | KS p=%.2e | RS p=%.2e', p_ks, p_rs));
legend(cohort_a_name, cohort_b_name);
box off;

subplot(2, 2, 4);
bar_x = [1 2];
bmn = [mean(d_er_a) mean(d_er_b)];
hold on;
bar(bar_x, bmn, 0.45, 'FaceColor', 'flat', 'CData', [col_a; col_b], 'EdgeColor', 'k');
errorbar(bar_x, bmn, [std(d_er_a)/sqrt(numel(d_er_a)), std(d_er_b)/sqrt(numel(d_er_b))], 'k.', 'LineWidth', 1.2);
set(gca, 'XTick', bar_x, 'XTickLabel', {cohort_a_name, cohort_b_name});
ylabel('Mean of per-neuron Δ event rate');
[p_er, ~] = ranksum(d_er_a, d_er_b);
title(sprintf('Pooled neurons: mean Δ event rate (± SE) | RS p=%.2e', p_er));
box off;

exportgraphics(gcf, fullfile(out_dir, 'fig_cohort_delta_population.png'), 'Resolution', 200);

figure('Color', 'w', 'Position', [30 30 1600 500]);
sgtitle({['Per-neuron Δ distributions (CDF) — pooled across mice.'], tw_label}, 'FontSize', 12, 'FontWeight', 'bold');
cdf_a = {d_mdf_a, d_er_a, d_std_a, d_sk_a};
cdf_b = {d_mdf_b, d_er_b, d_std_b, d_sk_b};
cdf_lbl = {'Δ Mean dF/F', 'Δ Event rate (1/s)', 'Δ Std dF/F', 'Δ Skewness'};
for i = 1:4
    subplot(1, 4, i);
    d1 = cdf_a{i};
    d2 = cdf_b{i};
    [f1, x1] = ecdf(d1);
    [f2, x2] = ecdf(d2);
    plot(x1, f1, 'Color', col_a, 'LineWidth', 2);
    hold on;
    plot(x2, f2, 'Color', col_b, 'LineWidth', 2);
    xlabel(cdf_lbl{i});
    ylabel('CDF');
    [~, p_ks] = kstest2(d1, d2);
    [p_rs, ~] = ranksum(d1, d2);
    title(sprintf('%s\nKS p=%.1e | RS p=%.1e', cdf_lbl{i}, p_ks, p_rs), 'FontSize', 9);
    xl = prctile([d1; d2], [0.5 99.5]);
    xlim(xl);
    legend(sprintf('%s med=%.3g', cohort_a_name, median(d1)), sprintf('%s med=%.3g', cohort_b_name, median(d2)), 'Location', 'southeast', 'FontSize', 7);
    box off;
end
exportgraphics(gcf, fullfile(out_dir, 'fig_cohort_delta_neuron_cdf.png'), 'Resolution', 200);

figure('Color', 'w', 'Position', [30 30 1600 650]);
sgtitle({sprintf('ΔZ heatmap: z_{inj}(t) - mean_t(z_{BL}(t)) | %d neurons sampled per cohort', max_show), tw_label}, ...
    'FontSize', 11, 'Interpreter', 'tex');
[Za, dur_a] = pooled_delta_z_matrix_legacy(Ba, max_show, 99);
[Zb, dur_b] = pooled_delta_z_matrix_legacy(Bb, max_show, 100);
c_hi = 5;
if ~isempty(Za) && ~isempty(Zb)
    c_hi = max([5; prctile(Za(:), 99.5); prctile(Zb(:), 99.5)]);
elseif ~isempty(Za)
    c_hi = max(5, prctile(Za(:), 99.5));
elseif ~isempty(Zb)
    c_hi = max(5, prctile(Zb(:), 99.5));
end
subplot(1, 2, 1);
if ~isempty(Za)
    imagesc([0 dur_a], [1 size(Za, 1)], Za);
    caxis([0 c_hi]);
    colormap(gca, cohort_cmap_inferno_legacy());
    xlabel('Time (s)');
    ylabel('Neuron (sorted by peak ΔZ)');
    title(sprintf('%s (pooled)', cohort_a_name));
end
subplot(1, 2, 2);
if ~isempty(Zb)
    imagesc([0 dur_b], [1 size(Zb, 1)], Zb);
    caxis([0 c_hi]);
    colormap(gca, cohort_cmap_inferno_legacy());
    cb = colorbar;
    cb.Label.String = 'ΔZ';
    xlabel('Time (s)');
    ylabel('Neuron (sorted by peak ΔZ)');
    title(sprintf('%s (pooled)', cohort_b_name));
end
exportgraphics(gcf, fullfile(out_dir, 'fig_cohort_delta_heatmaps.png'), 'Resolution', 200);

fp = params.frame_px;
[Ma_er, Ma_df, Mb_er, Mb_df] = mean_spatial_delta_maps_legacy(Ba, Bb, fp);
figure('Color', 'w', 'Position', [30 30 1500 650]);
sgtitle({'Mean spatial maps (cohort average over mice): injection-neuron Δ vs baseline session mean', tw_label}, 'FontSize', 10);

subplot(2, 2, 1);
if ~isempty(Ma_er)
    imagesc(Ma_er);
    axis image off;
    colormap(gca, cohort_cmap_diverging_legacy());
    caxis(cohort_symmetric_clim_legacy(Ma_er));
    colorbar;
    title(sprintf('%s: mean Δ event rate', cohort_a_name));
end

subplot(2, 2, 2);
if ~isempty(Mb_er)
    imagesc(Mb_er);
    axis image off;
    colormap(gca, cohort_cmap_diverging_legacy());
    caxis(cohort_symmetric_clim_legacy(Mb_er));
    colorbar;
    title(sprintf('%s: mean Δ event rate', cohort_b_name));
end

subplot(2, 2, 3);
if ~isempty(Ma_df)
    imagesc(Ma_df);
    axis image off;
    colormap(gca, cohort_cmap_diverging_legacy());
    caxis(cohort_symmetric_clim_legacy(Ma_df));
    colorbar;
    title(sprintf('%s: mean Δ mean dF/F', cohort_a_name));
end

subplot(2, 2, 4);
if ~isempty(Mb_df)
    imagesc(Mb_df);
    axis image off;
    colormap(gca, cohort_cmap_diverging_legacy());
    caxis(cohort_symmetric_clim_legacy(Mb_df));
    colorbar;
    title(sprintf('%s: mean Δ mean dF/F', cohort_b_name));
end

exportgraphics(gcf, fullfile(out_dir, 'fig_cohort_delta_spatial.png'), 'Resolution', 200);

fprintf('Saved figures + stats under: %s\n', out_dir);
end

function params = cohort_neuron_figure_params_legacy(time_window_mode)
params = struct( ...
    'z_thresh', 2, ...
    'time_bin', 1, ...
    'n_pc', 10, ...
    'n_corr_sample', 500, ...
    'frame_px', 2304, ...
    'skip_first_sec_bl', 0, ...
    'use_last_sec_bl', inf, ...
    'skip_first_sec_inj', 0, ...
    'use_last_sec_inj', inf);
switch lower(strtrim(time_window_mode))
    case 'matched'
        params.match_inj_last_to_baseline_dur = true;
    case 'whole'
        params.match_inj_last_to_baseline_dur = false;
    otherwise
        error('extract_cohort_neuron_delta_figures_legacy:BadTimeMode', ...
            'time_window_mode must be ''matched'' or ''whole'', got ''%s''.', time_window_mode);
end
end

function lim = cohort_symmetric_clim_legacy(M, default_abs)
if nargin < 2 || isempty(default_abs)
    default_abs = 1;
end
if isempty(M)
    lim = [-default_abs, default_abs];
    return;
end
v = double(M(:));
v = v(isfinite(v));
if isempty(v)
    lim = [-default_abs, default_abs];
    return;
end
sc = prctile(abs(v), 99);
if ~isfinite(sc) || sc < 1e-12
    sc = default_abs;
end
lim = [-sc, sc];
end

function [bundles, ids_ok] = collect_bundles_legacy(root, mouse_list, bl_f, inj_f, mirror, params)
bundles = {};
ids_ok = {};
for mi = 1:numel(mouse_list)
    mid = mouse_list{mi};
    bl_path  = fullfile(root, mid, bl_f, mirror);
    inj_path = fullfile(root, mid, inj_f, mirror);
    if exist(bl_path, 'dir') ~= 7 || exist(inj_path, 'dir') ~= 7
        warning('extract_cohort_neuron_delta_figures_legacy:Skip', 'Missing path for %s', mid);
        continue;
    end
    try
        B = extract_bl_inj_neuron_bundle(bl_path, inj_path, params);
        B.mouse_id = mid;
        B.mirror = mirror;
        bundles{end+1} = B; %#ok<AGROW>
        ids_ok{end+1} = mid; %#ok<AGROW>
        fprintf('  OK %s | n_bl=%d n_inj=%d\n', mid, B.n_neurons_bl, B.n_neurons_inj);
    catch ME
        warning('extract_cohort_neuron_delta_figures_legacy:BundleFail', '%s: %s', mid, ME.message);
    end
end
end

function log_cmp_legacy(fid, name, a, b)
[~, p_ks] = kstest2(a, b);
[p_rs, ~] = ranksum(a, b);
fprintf(fid, '%-40s  n_a=%6d n_b=%6d  KS p=%.4g  RS p=%.4g\n', name, numel(a), numel(b), p_ks, p_rs);
end

function [all_dpf_a, all_dpf_b, all_dpr_a, all_dpr_b, t_vec, mu_pf_a, sem_pf_a, mu_pf_b, sem_pf_b, mu_pr_a, sem_pr_a, mu_pr_b, sem_pr_b] = ...
    population_delta_series_legacy(Ba, Bb, time_bin)

all_dpf_a = [];
all_dpr_a = [];
for k = 1:numel(Ba)
    B = Ba{k};
    Tm = min(B.n_bins_bl, B.n_bins_inj);
    all_dpf_a = [all_dpf_a; B.pop_frac_inj(1:Tm) - B.pop_frac_bl(1:Tm)]; %#ok<AGROW>
    all_dpr_a = [all_dpr_a; B.pop_rate_inj(1:Tm) - B.pop_rate_bl(1:Tm)]; %#ok<AGROW>
end
all_dpf_b = [];
all_dpr_b = [];
for k = 1:numel(Bb)
    B = Bb{k};
    Tm = min(B.n_bins_bl, B.n_bins_inj);
    all_dpf_b = [all_dpf_b; B.pop_frac_inj(1:Tm) - B.pop_frac_bl(1:Tm)]; %#ok<AGROW>
    all_dpr_b = [all_dpr_b; B.pop_rate_inj(1:Tm) - B.pop_rate_bl(1:Tm)]; %#ok<AGROW>
end

Tm_a = inf;
for k = 1:numel(Ba)
    B = Ba{k};
    Tm_a = min(Tm_a, min(B.n_bins_bl, B.n_bins_inj));
end
Tm_b = inf;
for k = 1:numel(Bb)
    B = Bb{k};
    Tm_b = min(Tm_b, min(B.n_bins_bl, B.n_bins_inj));
end
Tm = min(Tm_a, Tm_b);
if ~isfinite(Tm) || Tm < 1
    t_vec = [];
    [mu_pf_a, sem_pf_a, mu_pf_b, sem_pf_b, mu_pr_a, sem_pr_a, mu_pr_b, sem_pr_b] = deal([]);
    return;
end

Matfa = nan(numel(Ba), Tm);
Matfb = nan(numel(Bb), Tm);
Matra = nan(numel(Ba), Tm);
Matrb = nan(numel(Bb), Tm);
for k = 1:numel(Ba)
    B = Ba{k};
    Matfa(k, :) = (B.pop_frac_inj(1:Tm) - B.pop_frac_bl(1:Tm))';
    Matra(k, :) = (B.pop_rate_inj(1:Tm) - B.pop_rate_bl(1:Tm))';
end
for k = 1:numel(Bb)
    B = Bb{k};
    Matfb(k, :) = (B.pop_frac_inj(1:Tm) - B.pop_frac_bl(1:Tm))';
    Matrb(k, :) = (B.pop_rate_inj(1:Tm) - B.pop_rate_bl(1:Tm))';
end
mu_pf_a = mean(Matfa, 1, 'omitnan');
sem_pf_a = std(Matfa, 0, 1, 'omitnan') ./ max(sqrt(sum(~isnan(Matfa), 1)), 1);
mu_pf_b = mean(Matfb, 1, 'omitnan');
sem_pf_b = std(Matfb, 0, 1, 'omitnan') ./ max(sqrt(sum(~isnan(Matfb), 1)), 1);
mu_pr_a = mean(Matra, 1, 'omitnan');
sem_pr_a = std(Matra, 0, 1, 'omitnan') ./ max(sqrt(sum(~isnan(Matra), 1)), 1);
mu_pr_b = mean(Matrb, 1, 'omitnan');
sem_pr_b = std(Matrb, 0, 1, 'omitnan') ./ max(sqrt(sum(~isnan(Matrb), 1)), 1);
t_vec = ((1:Tm)' - 0.5) * time_bin;
end

function [Zsort, dur_s] = pooled_delta_z_matrix_legacy(bundles, max_show, seed)
cols = [];
min_T = inf;
if nargin < 3
    seed = 99;
end
if isempty(bundles)
    Zsort = [];
    dur_s = 1;
    return;
end
for k = 1:numel(bundles)
    B = bundles{k};
    Tm = min(size(B.z_bl, 1), size(B.z_inj, 1));
    min_T = min(min_T, Tm);
end
if ~isfinite(min_T) || min_T < 2
    Zsort = [];
    dur_s = 1;
    return;
end
for k = 1:numel(bundles)
    B = bundles{k};
    zbl = B.z_bl(1:min_T, :);
    zij = B.z_inj(1:min_T, :);
    mu_t = mean(zbl, 2);
    dz = zij - mu_t;
    cols = [cols dz]; %#ok<AGROW>
end
n_take = min(max_show, size(cols, 2));
if n_take < 1
    Zsort = [];
    dur_s = min_T / bundles{1}.fs_inj;
    return;
end
rng(seed);
idx = randperm(size(cols, 2), n_take);
subm = cols(:, idx);
[~, pk] = max(subm, [], 1);
[~, ord] = sort(pk);
Zsort = subm(:, ord)';
dur_s = min_T / bundles{1}.fs_inj;
end

function [Ma_er, Ma_df, Mb_er, Mb_df] = mean_spatial_delta_maps_legacy(Ba, Bb, frame_px)
Ma_er = average_delta_map_legacy(Ba, frame_px, 'er');
Ma_df = average_delta_map_legacy(Ba, frame_px, 'mdf');
Mb_er = average_delta_map_legacy(Bb, frame_px, 'er');
Mb_df = average_delta_map_legacy(Bb, frame_px, 'mdf');
end

function M = average_delta_map_legacy(bundles, frame_px, which_m)
acc = zeros(frame_px, frame_px);
n = 0;
for k = 1:numel(bundles)
    B = bundles{k};
    switch which_m
        case 'er'
            v = B.event_rate_inj - mean(B.event_rate_bl);
        case 'mdf'
            v = B.mean_dff_inj - mean(B.mean_dff_bl);
        otherwise
            v = B.event_rate_inj - mean(B.event_rate_bl);
    end
    try
        m = project_inj_vector_to_map_legacy(B.S_inj_ir, B.S_inj_jc, B.S_inj_data, frame_px, v);
        acc = acc + m;
        n = n + 1;
    catch
    end
end
if n < 1
    M = [];
else
    M = acc / n;
end
end

function M = project_inj_vector_to_map_legacy(ir, jc, data, frame_px, vec)
n_pix = frame_px * frame_px;
n_neu = numel(vec);
col_idx = zeros(jc(end), 1);
for c = 1:n_neu
    col_idx(jc(c)+1:jc(c+1)) = c;
end
Ss = sparse(ir + 1, col_idx, data, n_pix, n_neu);
cs = full(sum(Ss, 1));
cs(cs == 0) = 1;
Sn = Ss ./ cs;
M = reshape(full(Sn * vec(:)), frame_px, frame_px);
end

function cmap = cohort_cmap_diverging_legacy()
half_n = 128;
cmap = [[linspace(0, 1, half_n) ones(1, half_n)]' ...
        [linspace(0, 1, half_n) linspace(1, 0, half_n)]' ...
        [ones(1, half_n) linspace(1, 0, half_n)]'];
end

function cmap = cohort_cmap_inferno_legacy()
t_cm = linspace(0, 1, 256)';
inf_r = min(1, max(0, -0.0155 + t_cm .* (5.27 + t_cm .* (-14.35 + t_cm .* (17.77 + t_cm .* (-9.42 + t_cm .* 1.76))))));
inf_g = min(1, max(0, -0.0123 + t_cm .* (-0.23 + t_cm .* (7.77 + t_cm .* (-17.05 + t_cm .* (16.95 + t_cm .* (-6.34)))))));
inf_b = min(1, max(0, 0.095 + t_cm .* (4.00 + t_cm .* (-15.24 + t_cm .* (23.56 + t_cm .* (-17.39 + t_cm .* 5.01))))));
cmap = [inf_r inf_g inf_b];
end
