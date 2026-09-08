function write_n18_clean_figs()
%WRITE_N18_CLEAN_FIGS  Cohort mean across 18 recordings (no day labels).

addpath(fileparts(mfilename("fullpath")));
root = fileparts(mfilename("fullpath"));
outDir = fullfile(root, "stim_triggered_figures", "combined_0813_0818", "n18");
if exist(outDir, "dir") ~= 7
    mkdir(outDir);
end

jobs = {
    fullfile(root, "stim_triggered_figures", "20250813_lab_meeting")
    fullfile(root, "stim_triggered_figures", "20250818_lab_meeting")
    };
results = {};
for j = 1:numel(jobs)
    T = readtable(fullfile(jobs{j}, "lab_meeting", "recording_index.csv"), "Delimiter", ",");
    for i = 1:height(T)
        movie = char(T.movie{i});
        tmp = load(fullfile(jobs{j}, movie, "stim_results.mat"), "R");
        results{end+1} = tmp.R; %#ok<AGROW>
    end
end

n = numel(results);
tq = -2:0.05:3;
psth = nan(n, numel(tq));
frac = nan(n, numel(tq));
staA = nan(n, numel(tq));
staE = nan(n, numel(tq));
staI = nan(n, numel(tq));
fe = nan(n, 1);
fi = nan(n, 1);
pp = nan(n, 1);
colExc = [0.85 0.20 0.15];
colInh = [0.10 0.35 0.85];
sem = @(x) std(x, 0, 1, "omitnan") ./ sqrt(max(sum(isfinite(x), 1), 1));

for i = 1:n
    R = results{i};
    psth(i, :) = interp1(R.t_peri, R.pop_psth, tq, "linear", NaN);
    frac(i, :) = interp1(R.t_peri, R.frac_psth, tq, "linear", NaN);
    if isfield(R, "sta_all"); staA(i, :) = interp1(R.t_peri, R.sta_all, tq, "linear", NaN); end
    if isfield(R, "sta_exc") && all(isfinite(R.sta_exc))
        staE(i, :) = interp1(R.t_peri, R.sta_exc, tq, "linear", NaN);
    end
    if isfield(R, "sta_inh") && all(isfinite(R.sta_inh))
        staI(i, :) = interp1(R.t_peri, R.sta_inh, tq, "linear", NaN);
    end
    fe(i) = R.frac_exc;
    fi(i) = R.frac_inh;
    pp(i) = R.p_pop;
end

mu = mean(psth, 1, "omitnan");
se = sem(psth);
muF = mean(frac, 1, "omitnan");
seF = sem(frac);

figure("Position", [20, 20, 1500, 560], "Color", "w");
sgtitle(sprintf("Population peri-poke response   (n = %d recordings)", n), ...
    "FontSize", 15, "FontWeight", "bold");
subplot(1, 2, 1);
plot(tq, psth', "Color", [0.78 0.78 0.78], "LineWidth", 0.8); hold on;
fill([tq fliplr(tq)], [mu-se fliplr(mu+se)], colExc, "FaceAlpha", 0.25, "EdgeColor", "none");
plot(tq, mu, "r", "LineWidth", 2.6);
xline(0, "k-", "LineWidth", 1.4);
xlabel("Time from poke (s)"); ylabel("Population mean z");
title("Gray = one recording   |   red = mean \pm SEM");
box off;
subplot(1, 2, 2);
plot(tq, frac', "Color", [0.78 0.78 0.78], "LineWidth", 0.8); hold on;
fill([tq fliplr(tq)], [muF-seF fliplr(muF+seF)], [0.3 0.3 0.3], "FaceAlpha", 0.2, "EdgeColor", "none");
plot(tq, muF, "k", "LineWidth", 2.6);
xline(0, "k-", "LineWidth", 1.4);
xlabel("Time from poke (s)"); ylabel("Fraction active (z>2)");
title("Population recruitment");
box off;
exportgraphics(gcf, fullfile(outDir, "fig1_n18_psth.png"), "Resolution", 200);

figure("Position", [20, 20, 1500, 560], "Color", "w");
sgtitle(sprintf("Poke-responsive cells   (n = %d recordings)", n), ...
    "FontSize", 15, "FontWeight", "bold");
subplot(1, 2, 1);
plot(tq, mean(staE, 1, "omitnan"), "Color", colExc, "LineWidth", 2.4); hold on;
plot(tq, mean(staI, 1, "omitnan"), "Color", colInh, "LineWidth", 2.4);
plot(tq, mean(staA, 1, "omitnan"), "k", "LineWidth", 2.0);
xline(0, "k-", "LineWidth", 1.4);
xlabel("Time from poke (s)"); ylabel("Mean STA (z)");
title("Excited vs inhibited vs all cells");
legend('Excited', 'Inhibited', 'All cells', 'Location', 'best');
box off;
subplot(1, 2, 2);
b = bar([1 2], [mean(fe) mean(fi)]*100, 0.55, "EdgeColor", "k"); hold on;
b.FaceColor = "flat"; b.CData = [colExc; colInh];
errorbar([1 2], [mean(fe) mean(fi)]*100, [std(fe) std(fi)]*100, "k.", "LineWidth", 1.4);
rng(0);
scatter(1 + 0.12*randn(n,1), 100*fe, 28, "k", "filled");
scatter(2 + 0.12*randn(n,1), 100*fi, 28, "k", "filled");
set(gca, "XTick", 1:2, "XTickLabel", {"Excited", "Inhibited"});
ylabel("% of neurons");
title(sprintf("%.1f%% excited,  %.1f%% inhibited   |   %d/%d recordings p<0.05", ...
    100*mean(fe), 100*mean(fi), sum(pp<0.05), n));
box off;
exportgraphics(gcf, fullfile(outDir, "fig2_n18_cells.png"), "Resolution", 200);
close all;
fprintf("n=18 figures -> %s\n", outDir);
end
