function extract_cohort_delta_pipeline_truncated()
%EXTRACT_COHORT_DELTA_PIPELINE_TRUNCATED  Same as extract_cohort_delta_pipeline, but metrics use:
%   Baseline = full trace after optional BL crops (default: no BL crop)
%   Injection = last N seconds where N = that baseline duration (match_inj_last_to_baseline_dur)
%
%   extract_cohort_delta_pipeline_truncated
%
% Cohort comparison defaults to PAIRED: cohort_a_mice{i} and cohort_b_mice{i} = same mouse.
% Statistics: paired t-test and Wilcoxon signed-rank on Δ (inj−baseline) per cohort.
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
mirror_subfolders = {'mirror1path', 'mirror2path'};

out_dir = 'C:\Users\hsollim\Documents\MATLAB\MATLAB\2p\cohort_delta_melittin_vs_morphine_INJlastMatchBL';

paired_cohort_compare = true;  % row i in A ↔ row i in B (same animal)

params = struct( ...
    'z_thresh', 2, ...
    'time_bin', 1, ...
    'n_pc', 10, ...
    'n_corr_sample', 500, ...
    'skip_first_sec_bl', 0, ...
    'use_last_sec_bl', inf, ...
    'skip_first_sec_inj', 0, ...
    'use_last_sec_inj', inf, ...
    'match_inj_last_to_baseline_dur', true);
%% ========================================================================

close all;
clc;

if exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end

mirror_list = extract_cohort_normalize_mirror_list(mirror_subfolders);

if paired_cohort_compare
    rows = collect_rows_paired_trunc(cohort_a_name, cohort_a_root, cohort_a_mice, ...
        cohort_b_name, cohort_b_root, cohort_b_mice, ...
        baseline_folder, injection_folder, mirror_list, params);
else
    rows = [ ...
        collect_rows_one_cohort_trunc(1, cohort_a_name, cohort_a_root, cohort_a_mice, baseline_folder, injection_folder, mirror_list, params); ...
        collect_rows_one_cohort_trunc(2, cohort_b_name, cohort_b_root, cohort_b_mice, baseline_folder, injection_folder, mirror_list, params)];
end

if isempty(rows)
    error('No mice processed; check paths and file resolution.');
end

metric_names = get_delta_metric_names_trunc();
Rw = struct2table(rows);
writetable(Rw, fullfile(out_dir, 'cohort_delta_per_mouse.csv'));

g1 = Rw.cohort == 1;
g2 = Rw.cohort == 2;
fid = fopen(fullfile(out_dir, 'delta_group_comparison.txt'), 'w');
fprintf(fid, 'Delta = injection − baseline | Baseline: full (after BL time params) | Injection: LAST N s where N = baseline duration per mouse\n');
fprintf(fid, 'Cohort 1: %s | Cohort 2: %s\n', cohort_a_name, cohort_b_name);
fprintf(fid, 'Focal planes (mirrors): %s\n', strjoin(mirror_list, ', '));
fprintf(fid, 'Figure fig_cohort_delta_comparison.png: when paired, one row per mouse (pair_idx); columns = metrics.\n');
if paired_cohort_compare
    fprintf(fid, 'Comparison: PAIRED — same mouse row i and same mirror (paired t-test + signrank).\n');
else
    fprintf(fid, 'Comparison: INDEPENDENT cohorts. Wilcoxon rank-sum.\n');
end
fprintf(fid, 'Note: small n at mouse×plane level; exploratory.\n\n');

p_fields = metric_names;
p_paired_t = nan(size(p_fields));
p_signrank = nan(size(p_fields));
p_ranksum = nan(size(p_fields));

for k = 1:numel(p_fields)
    f = p_fields{k};
    if ~ismember(f, Rw.Properties.VariableNames)
        continue;
    end
    if paired_cohort_compare
        T1 = Rw(Rw.cohort == 1, :);
        T2 = Rw(Rw.cohort == 2, :);
        [T1, T2] = extract_cohort_align_paired_tables(T1, T2);
        x1 = T1.(f);
        x2 = T2.(f);
        [p_paired_t(k), p_signrank(k), p_ranksum(k)] = extract_cohort_delta_pair_stats(x1, x2, true);
        fprintf(fid, '%-28s  d %s = [%s]  d %s = [%s]  paired_t p=%.4g  signrank p=%.4g\n', ...
            f, cohort_a_name, num2str(x1', '%.4g '), cohort_b_name, num2str(x2', '%.4g '), ...
            p_paired_t(k), p_signrank(k));
    else
        x1 = Rw.(f)(g1);
        x2 = Rw.(f)(g2);
        [p_paired_t(k), p_signrank(k), p_ranksum(k)] = extract_cohort_delta_pair_stats(x1, x2, false);
        fprintf(fid, '%-28s  d cohort1 = [%s]  d cohort2 = [%s]  ranksum p = %.4g\n', ...
            f, num2str(x1', '%.4g '), num2str(x2', '%.4g '), p_ranksum(k));
    end
end
fclose(fid);

if paired_cohort_compare
    sgtitle_str = sprintf('\\Delta (inj−baseline) | PAIRED | INJ last = BL duration | %s vs %s', ...
        cohort_a_name, cohort_b_name);
    title_fs = 10;
else
    sgtitle_str = sprintf('\\Delta (inj−baseline) | INJ last = BL duration | %s vs %s', ...
        cohort_a_name, cohort_b_name);
    title_fs = 11;
end
extract_cohort_plot_fig_delta_comparison(out_dir, Rw, paired_cohort_compare, cohort_a_name, cohort_b_name, ...
    p_fields, p_paired_t, p_signrank, p_ranksum, sgtitle_str, @cohort_dot_label_trunc, title_fs);

fprintf('\nDone (truncated windows). Outputs in:\n  %s\n', out_dir);
fprintf('Elapsed %.1f s\n', toc);

end

function rows = collect_rows_paired_trunc(name_a, root_a, mice_a, name_b, root_b, mice_b, bl_f, inj_f, mirror_list, params)
rows = [];
n = numel(mice_a);
if numel(mice_b) ~= n
    error('extract_cohort_delta_pipeline_truncated:PairedMouseCount', ...
        'paired_cohort_compare requires equal mouse list lengths; got %d vs %d.', n, numel(mice_b));
end
for i = 1:n
    ma = mice_a{i};
    mb = mice_b{i};
    for km = 1:numel(mirror_list)
        mirror = mirror_list{km};
        bl_a  = fullfile(root_a, ma, bl_f, mirror);
        inj_a = fullfile(root_a, ma, inj_f, mirror);
        bl_b  = fullfile(root_b, mb, bl_f, mirror);
        inj_b = fullfile(root_b, mb, inj_f, mirror);
        if exist(bl_a, 'dir') ~= 7 || exist(inj_a, 'dir') ~= 7 || exist(bl_b, 'dir') ~= 7 || exist(inj_b, 'dir') ~= 7
            warning('extract_cohort_delta_pipeline_truncated:SkipPair', ...
                'Skip pair %d mirror %s (%s vs %s).', i, mirror, ma, mb);
            continue;
        end
        fprintf('--- pair %d | %s | %s: %s  ↔  %s: %s ---\n', i, mirror, name_a, ma, name_b, mb);
        try
            D_a = extract_cohort_delta_from_session_pair(bl_a, inj_a, params);
            D_b = extract_cohort_delta_from_session_pair(bl_b, inj_b, params);
        catch ME
            warning('extract_cohort_delta_pipeline_truncated:PairFailed', 'Pair %d mirror %s: %s', i, mirror, ME.message);
            continue;
        end
        rA = struct('pair_idx', i, 'mirror', {mirror}, 'cohort', 1, 'cohort_name', {name_a}, 'mouse_id', {ma});
        rA = merge_row_trunc(rA, D_a);
        rB = struct('pair_idx', i, 'mirror', {mirror}, 'cohort', 2, 'cohort_name', {name_b}, 'mouse_id', {mb});
        rB = merge_row_trunc(rB, D_b);
        rows = [rows; rA; rB]; %#ok<AGROW>
    end
end
end

function rows = collect_rows_one_cohort_trunc(cohort_id, cohort_name, root, mouse_list, bl_f, inj_f, mirror_list, params)
rows = [];
for mi = 1:numel(mouse_list)
    mouse_id = mouse_list{mi};
    for km = 1:numel(mirror_list)
        mirror = mirror_list{km};
        bl_path  = fullfile(root, mouse_id, bl_f, mirror);
        inj_path = fullfile(root, mouse_id, inj_f, mirror);
        if exist(bl_path, 'dir') ~= 7 || exist(inj_path, 'dir') ~= 7
            warning('extract_cohort_delta_pipeline_truncated:SkipMouse', 'Missing %s / %s / %s', mouse_id, cohort_name, mirror);
            continue;
        end
        fprintf('--- %s | %s | %s ---\n', cohort_name, mouse_id, mirror);
        try
            D = extract_cohort_delta_from_session_pair(bl_path, inj_path, params);
        catch ME
            warning('extract_cohort_delta_pipeline_truncated:MetricsFailed', '%s %s: %s', mouse_id, mirror, ME.message);
            continue;
        end
        r = struct('pair_idx', nan, 'mirror', {mirror}, 'cohort', cohort_id, 'cohort_name', {cohort_name}, 'mouse_id', {mouse_id});
        r = merge_row_trunc(r, D);
        rows = [rows; r]; %#ok<AGROW>
    end
end
end

function lbl = cohort_dot_label_trunc(T, pi)
mid = char(string(T.mouse_id(pi)));
if ismember('mirror', T.Properties.VariableNames)
    lbl = sprintf('%s | %s', mid, char(string(T.mirror(pi))));
else
    lbl = mid;
end
end

function m = merge_row_trunc(a, b)
fn = fieldnames(b);
m = a;
for kk = 1:numel(fn)
    m.(fn{kk}) = b.(fn{kk});
end
end

function fn = get_delta_metric_names_trunc()
fn = {'d_mean_mean_dff', 'd_median_mean_dff', 'd_mean_std_dff', 'd_median_std_dff', ...
    'd_mean_skew_dff', 'd_median_skew_dff', 'd_mean_event_rate', 'd_median_event_rate', ...
    'd_mean_pop_frac', 'd_mean_pop_rate_z', 'd_median_pw_corr', 'd_pc1_var_pct', 'd_pc1_to_pc3_var_pct', ...
    'd_n_neurons'};
end
