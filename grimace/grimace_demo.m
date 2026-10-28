function grimace_demo()
%GRIMACE_DEMO  Synthetic PainFace-like CSV + figures (no video required).

here = fileparts(mfilename("fullpath"));
addpath(here);
outDir = fullfile(here, "demo_output");
if exist(outDir, "dir") ~= 7
    mkdir(outDir);
end

rng(1);
t = (0:1:599)';
n = numel(t);
pokes = 8:12:580;
mgs = 1.2 + 0.8*randn(n, 1);
for k = 1:numel(pokes)
    d = t - pokes(k);
    bump = 3.5 * exp(-(d - 0.6).^2 / (2*0.35^2));
    mgs = mgs + bump;
end
mgs = min(max(round(mgs), 0), 8);
nFau = 4 * ones(n, 1);
nFau(1:20:end) = 2;

csvFile = fullfile(outDir, "example_painface.csv");
writetable(table(t, (1:n)', mgs, nFau, VariableNames=["time_s", "frame", "mgs", "n_fau"]), ...
    csvFile);

G = grimace_load(csvFile, "Source", "painface");
M = grimace_metrics(G, "MinMaxFau", 50);
S = grimace_align_pokes(G, pokes);
grimace_plot(G, M, S, outDir);
fprintf("demo mean MGS %.2f, high-grimace %.1f%%, peri Delta %.2f p=%.3g\n", ...
    M.mean_mgs, M.pct_high, S.mean_amp, S.p_amp);
end
