function extract_cohort_delta_pipeline_truncated()
%EXTRACT_COHORT_DELTA_PIPELINE_TRUNCATED  Same as extract_cohort_delta_pipeline, but metrics use:
%   Baseline = full recording (no crop)
%   Injection = last 600 s only (same logic as extract_analysis2_truncated)
%
%   extract_cohort_delta_pipeline_truncated
%
% Outputs go to a different folder than the "whole recording" pipeline so results are not overwritten.

tic;

%% ========================= CONFIG =======================================
cohort_a_name = 'melittin_only';
cohort_a_root = '\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon';
cohort_a_mice = {'CS1014-1b', 'CS1014-2w'};

cohort_b_name = 'melittin_morphine';
cohort_b_root = '\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\morphine_20260408';
cohort_b_mice = {'1014m01b', '1014m02w'};

baseline_folder = 'Baseline';
injection_folder = 'Injection';
mirror_subfolder = 'mirror1path';

inj_last_sec = 600;

out_dir = 'C:\Users\hsollim\Documents\MATLAB\MATLAB\2p\cohort_delta_melittin_vs_morphine_BLfull_INJlast600';

params = struct( ...
    'z_thresh', 2, ...
    'time_bin', 1, ...
    'n_pc', 10, ...
    'n_corr_sample', 500, ...
    'skip_first_sec_bl', 0, ...
    'use_last_sec_bl', inf, ...
    'skip_first_sec_inj', 0, ...
    'use_last_sec_inj', inj_last_sec);
%% ========================================================================

close all;
clc;

if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

rows = [];

rows = [rows; process_cohort_trunc(1, cohort_a_name, cohort_a_root, cohort_a_mice, baseline_folder, injection_folder, mirror_subfolder, params)]; %#ok<AGROW>
rows = [rows; process_cohort_trunc(2, cohort_b_name, cohort_b_root, cohort_b_mice, baseline_folder, injection_folder, mirror_subfolder, params)];

if isempty(rows)
    error('No mice processed; check paths and file resolution.');
end

metric_names = get_delta_metric_names_trunc();
Rw = struct2table(rows);
writetable(Rw, fullfile(out_dir, 'cohort_delta_per_mouse.csv'));

g1 = Rw.cohort == 1;
g2 = Rw.cohort == 2;
fid = fopen(fullfile(out_dir, 'delta_group_comparison.txt'), 'w');
fprintf(fid, 'Delta = injection − baseline | Baseline FULL | Injection LAST %d s only\n', inj_last_sec);
fprintf(fid, 'Cohort 1: %s | Cohort 2: %s\n', cohort_a_name, cohort_b_name);
fprintf(fid, 'Mirror: %s\n', mirror_subfolder);
fprintf(fid, 'Note: ranksum with n=2 vs 2 mice is extremely low power; use for exploration only.\n\n');

p_fields = metric_names;
p_vals = nan(size(p_fields));
for k = 1:numel(p_fields)
    f = p_fields{k};
    if ~ismember(f, Rw.Properties.VariableNames)
        continue;
    end
    x1 = Rw.(f)(g1);
    x2 = Rw.(f)(g2);
    if numel(x1) >= 1 && numel(x2) >= 1
        [p_vals(k), ~] = ranksum(x1, x2);
    end
    fprintf(fid, '%-28s  d cohort1 = [%s]  d cohort2 = [%s]  ranksum p = %.4g\n', ...
        f, num2str(x1', '%.4g '), num2str(x2', '%.4g '), p_vals(k));
end
fclose(fid);

plot_metrics = {'d_mean_event_rate', 'd_mean_mean_dff', 'd_mean_pop_frac', 'd_median_pw_corr', 'd_pc1_var_pct', 'd_pc1_to_pc3_var_pct'};
plot_metrics = plot_metrics(ismember(plot_metrics, Rw.Properties.VariableNames));
nP = numel(plot_metrics);
if nP > 0
    figure('Color', 'w', 'Position', [40 40 200 * nP 400]);
    sgtitle(sprintf('\\Delta (inj−baseline) | BL full, INJ last %ds | %s vs %s', inj_last_sec, cohort_a_name, cohort_b_name), 'FontSize', 11);
    col1 = [0.13 0.59 0.95];
    col2 = [1.00 0.34 0.13];
    for k = 1:nP
        subplot(1, nP, k);
        hold on;
        y1 = Rw.(plot_metrics{k})(g1);
        y2 = Rw.(plot_metrics{k})(g2);
        x1 = 1 + 0.08 * (rand(sum(g1), 1) - 0.5);
        x2 = 2 + 0.08 * (rand(sum(g2), 1) - 0.5);
        scatter(x1, y1, 64, col1, 'filled', 'MarkerEdgeColor', 'k');
        scatter(x2, y2, 64, col2, 'filled', 'MarkerEdgeColor', 'k');
        yy = [y1(:); y2(:)];
        if ~isempty(yy)
            plot([1 1], [min(yy) max(yy)], 'k:', 'LineWidth', 0.5);
            plot([2 2], [min(yy) max(yy)], 'k:', 'LineWidth', 0.5);
        end
        set(gca, 'XTick', [1 2], 'XTickLabel', {cohort_a_name, cohort_b_name}, 'XTickLabelRotation', 25);
        ylabel(strrep(plot_metrics{k}, '_', '\_'));
        pk = find(strcmp(p_fields, plot_metrics{k}), 1);
        if ~isempty(pk)
            title(sprintf('p=%.3g', p_vals(pk)));
        end
        box off;
    end
    exportgraphics(gcf, fullfile(out_dir, 'fig_cohort_delta_comparison.png'), 'Resolution', 200);
end

fprintf('\nDone (truncated INJ window). Outputs in:\n  %s\n', out_dir);
fprintf('Elapsed %.1f s\n', toc);

end

function rows = process_cohort_trunc(cohort_id, cohort_name, root, mouse_list, bl_f, inj_f, mirror, params)
rows = [];
for mi = 1:numel(mouse_list)
    mouse_id = mouse_list{mi};
    bl_path  = fullfile(root, mouse_id, bl_f, mirror);
    inj_path = fullfile(root, mouse_id, inj_f, mirror);
    if exist(bl_path, 'dir') ~= 7 || exist(inj_path, 'dir') ~= 7
        warning('extract_cohort_delta_pipeline_truncated:SkipMouse', 'Missing folder for %s / %s', mouse_id, cohort_name);
        continue;
    end
    fprintf('--- %s | %s ---\n', cohort_name, mouse_id);
    try
        m_bl = extract_session_metrics(bl_path, params, 'baseline');
        m_inj = extract_session_metrics(inj_path, params, 'injection');
    catch ME
        warning('extract_cohort_delta_pipeline_truncated:MetricsFailed', '%s: %s', mouse_id, ME.message);
        continue;
    end
    D = compute_delta_metrics_trunc(m_inj, m_bl);
    r = struct('cohort', cohort_id, 'cohort_name', {cohort_name}, 'mouse_id', {mouse_id});
    r = merge_struct_trunc(r, D);
    rows = [rows; r]; %#ok<AGROW>
end
end

function D = compute_delta_metrics_trunc(m_inj, m_bl)
fn = {'mean_mean_dff', 'median_mean_dff', 'mean_std_dff', 'median_std_dff', ...
    'mean_skew_dff', 'median_skew_dff', 'mean_event_rate', 'median_event_rate', ...
    'mean_pop_frac', 'mean_pop_rate_z', 'median_pw_corr', 'pc1_var_pct', 'pc1_to_pc3_var_pct'};
D = struct();
for k = 1:numel(fn)
    f = fn{k};
    D.(['d_' f]) = m_inj.(f) - m_bl.(f);
end
D.d_n_neurons = m_inj.n_neurons - m_bl.n_neurons;
D.bl_n_neurons = m_bl.n_neurons;
D.inj_n_neurons = m_inj.n_neurons;
end

function m = merge_struct_trunc(a, b)
fn = fieldnames(b);
m = a;
for k = 1:numel(fn)
    m.(fn{k}) = b.(fn{k});
end
end

function fn = get_delta_metric_names_trunc()
fn = {'d_mean_mean_dff', 'd_median_mean_dff', 'd_mean_std_dff', 'd_median_std_dff', ...
    'd_mean_skew_dff', 'd_median_skew_dff', 'd_mean_event_rate', 'd_median_event_rate', ...
    'd_mean_pop_frac', 'd_mean_pop_rate_z', 'd_median_pw_corr', 'd_pc1_var_pct', 'd_pc1_to_pc3_var_pct', ...
    'd_n_neurons'};
end
