function extract_cohort_plot_fig_delta_comparison(out_dir, Rw, paired_cohort_compare, cohort_a_name, cohort_b_name, ...
    p_fields, p_paired_t, p_signrank, p_ranksum, sgtitle_str, dot_label_fun, title_font_size)
%EXTRACT_COHORT_PLOT_FIG_DELTA_COMPARISON  Scatter + paired lines per metric; optional row per mouse (pair_idx).
%
%   When paired_cohort_compare and table has finite pair_idx, layout is nMouse × nP (one row per mouse,
%   columns = metrics). Global p-values (full n) shown in each subplot title. Y-limits matched per column.

if nargin < 12 || isempty(title_font_size)
    title_font_size = 11;
end

plot_metrics = {'d_mean_event_rate', 'd_mean_mean_dff', 'd_mean_pop_frac', 'd_median_pw_corr', 'd_pc1_var_pct', 'd_pc1_to_pc3_var_pct'};
plot_metrics = plot_metrics(ismember(plot_metrics, Rw.Properties.VariableNames));
nP = numel(plot_metrics);
if nP == 0
    return;
end

T1plot = Rw(Rw.cohort == 1, :);
T2plot = Rw(Rw.cohort == 2, :);
if paired_cohort_compare
    [T1plot, T2plot] = extract_cohort_align_paired_tables(T1plot, T2plot);
end

split_by_mouse = paired_cohort_compare && ismember('pair_idx', T1plot.Properties.VariableNames) ...
    && any(isfinite(T1plot.pair_idx(:)));
if split_by_mouse
    mouse_ids = unique(T1plot.pair_idx(isfinite(T1plot.pair_idx)));
    mouse_ids = sort(mouse_ids(:)');
    nMouse = numel(mouse_ids);
else
    nMouse = 1;
    mouse_ids = [];
end

figH = 430;
if split_by_mouse && nMouse > 1
    figH = min(220 + 270 * nMouse, 1200);
end
figure('Color', 'w', 'Position', [40 40 min(240 * nP, 2400) figH]);
sgtitle(sgtitle_str, 'FontSize', title_font_size);

col1 = [0.13 0.59 0.95];
col2 = [1.00 0.34 0.13];
rng(1);

ylim_lo = inf(1, nP);
ylim_hi = -inf(1, nP);

for mi = 1:nMouse
    if split_by_mouse && ~isempty(mouse_ids)
        mid = mouse_ids(mi);
        idx = T1plot.pair_idx == mid & isfinite(T1plot.pair_idx);
        T1 = T1plot(idx, :);
        T2 = T2plot(idx, :);
        if size(T1, 1) > 0 && ismember('mouse_id', T1.Properties.VariableNames)
            row_lbl = sprintf('%s ↔ %s', char(string(T1.mouse_id(1))), char(string(T2.mouse_id(1))));
        elseif size(T1, 1) > 0
            row_lbl = sprintf('Pair %g', mid);
        else
            row_lbl = sprintf('Pair %g (no rows)', mid);
        end
    else
        T1 = T1plot;
        T2 = T2plot;
        row_lbl = '';
    end

    for kk = 1:nP
        subplot(nMouse, nP, (mi - 1) * nP + kk);
        hold on;
        f = plot_metrics{kk};
        y1 = T1.(f);
        y2 = T2.(f);
        xj1 = 1 + 0.06 * (rand(numel(y1), 1) - 0.5);
        xj2 = 2 + 0.06 * (rand(numel(y2), 1) - 0.5);
        if paired_cohort_compare && numel(y1) == numel(y2) && numel(y1) > 0
            for pi = 1:numel(y1)
                plot([xj1(pi) xj2(pi)], [y1(pi) y2(pi)], 'Color', [0.65 0.65 0.65], 'LineWidth', 0.9);
            end
        end
        scatter(xj1, y1, 64, col1, 'filled', 'MarkerEdgeColor', 'k');
        scatter(xj2, y2, 64, col2, 'filled', 'MarkerEdgeColor', 'k');
        for pi = 1:numel(y1)
            text(xj1(pi) + 0.05, y1(pi), dot_label_fun(T1, pi), ...
                'FontSize', 6, 'Interpreter', 'none', 'VerticalAlignment', 'middle');
        end
        for pi = 1:numel(y2)
            text(xj2(pi) - 0.05, y2(pi), dot_label_fun(T2, pi), ...
                'FontSize', 6, 'Interpreter', 'none', 'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'middle');
        end
        yy = [y1(:); y2(:)];
        if ~isempty(yy)
            plot([1 1], [min(yy) max(yy)], 'k:', 'LineWidth', 0.5);
            plot([2 2], [min(yy) max(yy)], 'k:', 'LineWidth', 0.5);
        end
        set(gca, 'XTick', [1 2], 'XTickLabel', {cohort_a_name, cohort_b_name}, 'XTickLabelRotation', 25);
        if split_by_mouse && ~isempty(strtrim(row_lbl)) && kk == 1
            ylabel(sprintf('%s\n%s', row_lbl, strrep(f, '_', '\_')));
        else
            ylabel(strrep(f, '_', '\_'));
        end
        pk = find(strcmp(p_fields, f), 1);
        if ~isempty(pk)
            if paired_cohort_compare
                title(sprintf('p_t=%.3g | p_{sr}=%.3g', p_paired_t(pk), p_signrank(pk)), 'FontSize', 9);
            else
                title(sprintf('p=%.3g', p_ranksum(pk)), 'FontSize', 9);
            end
        end
        box off;
        if ~isempty(yy)
            yl = ylim;
            ylim_lo(kk) = min(ylim_lo(kk), yl(1));
            ylim_hi(kk) = max(ylim_hi(kk), yl(2));
        end
    end
end

for mi = 1:nMouse
    for kk = 1:nP
        if isfinite(ylim_lo(kk)) && isfinite(ylim_hi(kk)) && ylim_hi(kk) > ylim_lo(kk)
            subplot(nMouse, nP, (mi - 1) * nP + kk);
            ylim([ylim_lo(kk), ylim_hi(kk)]);
        end
    end
end

exportgraphics(gcf, fullfile(out_dir, 'fig_cohort_delta_comparison.png'), 'Resolution', 200);
end
