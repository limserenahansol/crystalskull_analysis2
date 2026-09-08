function [p_paired_t, p_signrank, p_ranksum] = extract_cohort_delta_pair_stats(x_a, x_b, paired)
%EXTRACT_COHORT_DELTA_PAIR_STATS  Compare cohort A vs B on a metric (vectors of Δ per mouse).
%
%   paired — if true, x_a(i) and x_b(i) are the same biological mouse (row-aligned).
%            Paired t-test (ttest) and Wilcoxon signed-rank (signrank) on paired differences.
%            Requires numel(x_a)==numel(x_b). For n<2, paired t p is NaN.
%   paired — if false, independent samples: Wilcoxon rank-sum (ranksum).

x_a = x_a(:);
x_b = x_b(:);
p_paired_t = nan;
p_signrank = nan;
p_ranksum = nan;

if paired
    if numel(x_a) ~= numel(x_b)
        error('extract_cohort_delta_pair_stats:PairedLengthMismatch', ...
            'Paired mode requires equal n; got %d vs %d.', numel(x_a), numel(x_b));
    end
    n = numel(x_a);
    if n >= 2
        try
            [~, p_paired_t] = ttest(x_a, x_b);
        catch
            p_paired_t = nan;
        end
    end
    if n >= 1
        try
            p_signrank = signrank(x_a, x_b);
        catch
            p_signrank = nan;
        end
    end
else
    if numel(x_a) >= 1 && numel(x_b) >= 1
        try
            p_ranksum = ranksum(x_a, x_b);
        catch
            p_ranksum = nan;
        end
    end
end

end
