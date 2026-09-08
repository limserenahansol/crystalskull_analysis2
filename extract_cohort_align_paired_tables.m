function [T1o, T2o] = extract_cohort_align_paired_tables(T1, T2)
%EXTRACT_COHORT_ALIGN_PAIRED_TABLES  Align cohort 1 vs 2 rows for paired tests/plots.
%   Match on (pair_idx, mirror). Rows present in only one cohort are dropped.
%   If variable mirror is absent, match on pair_idx only (legacy single-plane).

if isempty(T1) || isempty(T2)
    T1o = T1;
    T2o = T2;
    return;
end

if ismember('mirror', T1.Properties.VariableNames) && ismember('mirror', T2.Properties.VariableNames)
    T1s = sortrows(T1, {'pair_idx', 'mirror'});
    T2s = sortrows(T2, {'pair_idx', 'mirror'});
    k1 = strcat(string(T1s.pair_idx), '|', string(T1s.mirror));
    k2 = strcat(string(T2s.pair_idx), '|', string(T2s.mirror));
else
    T1s = sortrows(T1, 'pair_idx');
    T2s = sortrows(T2, 'pair_idx');
    k1 = string(T1s.pair_idx);
    k2 = string(T2s.pair_idx);
end

k1c = cellstr(k1);
k2c = cellstr(k2);
keys = intersect(k1c, k2c);
keys = sort(keys(:)');
n = numel(keys);
if n == 0
    T1o = T1([], :);
    T2o = T2([], :);
    return;
end
idx1 = zeros(n, 1);
idx2 = zeros(n, 1);
for i = 1:n
    ki = keys{i};
    idx1(i) = find(strcmp(k1c, ki), 1);
    idx2(i) = find(strcmp(k2c, ki), 1);
end
T1o = T1s(idx1, :);
T2o = T2s(idx2, :);

end
