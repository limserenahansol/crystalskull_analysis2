function extract_stim_pool_days()
%EXTRACT_STIM_POOL_DAYS  Average peri-poke traces across 13 Aug + 18 Aug.
%
% Loads saved stim_results.mat (no re-EXTRACT). Each movie is one trace;
% the thick black line is the mean of all 18 movies. Day means are overlaid
% because t=0 is a brightness peak on 13 Aug and an intensity dip on 18 Aug.

addpath(fileparts(mfilename("fullpath")));
root = fileparts(mfilename("fullpath"));
outRoot = fullfile(root, "stim_triggered_figures", "combined_0813_0818");
avgDir = fullfile(outRoot, "pooled_average");
if exist(avgDir, "dir") ~= 7
    mkdir(avgDir);
end

jobs = {
    fullfile(root, "stim_triggered_figures", "20250813_lab_meeting"), "0813", "peak"
    fullfile(root, "stim_triggered_figures", "20250818_lab_meeting"), "0818", "dip"
    };

results = {};
days = {};
labels = {};
alignFiles = {};
for j = 1:size(jobs, 1)
    peakRoot = jobs{j, 1};
    day = jobs{j, 2};
    idxFile = fullfile(peakRoot, "lab_meeting", "recording_index.csv");
    T = readtable(idxFile, "Delimiter", ",");
    for i = 1:height(T)
        movie = char(T.movie{i});
        label = char(T.label{i});
        resFile = fullfile(peakRoot, movie, "stim_results.mat");
        if exist(resFile, "file") ~= 2
            warning("extract_stim_pool_days:Missing", "No stim_results.mat for %s", movie);
            continue
        end
        tmp = load(resFile, "R");
        R = tmp.R;
        R.movie = movie;
        R.day = day;
        R.label = label;
        results{end+1} = R; %#ok<AGROW>
        days{end+1} = day; %#ok<AGROW>
        labels{end+1} = sprintf("%s_%s", day, label); %#ok<AGROW>
        af = fullfile(peakRoot, "lab_meeting", "aligned_per_recording", ...
            sprintf("%s_fig4_aligned.png", label));
        if exist(af, "file") ~= 2
            af = fullfile(peakRoot, "lab_meeting", "aligned_per_recording", ...
                "cohort_aligned_sig.png");
        end
        alignFiles{end+1} = af; %#ok<AGROW>
        fprintf("  loaded %s  %s  n=%d  pokes=%d\n", day, label, R.n_neurons, R.n_pokes_total);
    end
end

params = struct("peri_win", [-2, 3], "z_thresh", 2);
write_pooled_figures(results, days, labels, alignFiles, avgDir, params);
write_pooled_csv(results, days, labels, outRoot);
fprintf("Pooled pack -> %s\n", avgDir);
end

function write_pooled_csv(results, days, labels, outRoot)
n = numel(results);
csvPath = fullfile(outRoot, "pooled_stim_summary.csv");
fid = fopen(csvPath, "w");
fprintf(fid, ['day,label,movie,event_field,n_neurons,n_pokes,n_trials,' ...
    'n_exc,frac_exc,n_inh,frac_inh,mean_amp,p_pop,iti_median\n']);
for i = 1:n
    R = results{i};
    fprintf(fid, "%s,%s,%s,%s,%d,%d,%d,%d,%.4f,%d,%.4f,%.5f,%.3e,%.3f\n", ...
        days{i}, labels{i}, R.movie, R.event_field, R.n_neurons, R.n_pokes_total, ...
        R.n_trials, R.n_exc, R.frac_exc, R.n_inh, R.frac_inh, R.mean_amp, ...
        R.p_pop, R.iti_median);
end
fclose(fid);
end

function write_pooled_figures(results, days, labels, alignFiles, avgDir, params)
n = numel(results);
is13 = false(1, n);
is18 = false(1, n);
for ii = 1:n
    d = char(days{ii});
    is13(ii) = strcmp(d, "0813") || strcmp(d, '0813');
    is18(ii) = strcmp(d, "0818") || strcmp(d, '0818');
end
fprintf("Pooled movies: %d  (0813=%d, 0818=%d)\n", n, sum(is13), sum(is18));
tq = params.peri_win(1):0.05:params.peri_win(2);
psth = nan(n, numel(tq));
frac = nan(n, numel(tq));
staA = nan(n, numel(tq));
staE = nan(n, numel(tq));
staI = nan(n, numel(tq));
wN = nan(n, 1);
wT = nan(n, 1);
fe = nan(n, 1);
fi = nan(n, 1);
pp = nan(n, 1);
npk = nan(n, 1);
nn = nan(n, 1);
nexc = nan(n, 1);
ninh = nan(n, 1);
alli = [];
col13 = [0.85 0.45 0.10];
col18 = [0.10 0.55 0.62];
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
    wN(i) = R.n_neurons;
    wT(i) = R.n_trials;
    fe(i) = R.frac_exc;
    fi(i) = R.frac_inh;
    pp(i) = R.p_pop;
    npk(i) = R.n_pokes_total;
    nn(i) = R.n_neurons;
    nexc(i) = R.n_exc;
    ninh(i) = R.n_inh;
    if isfield(R, "iti") && ~isempty(R.iti)
        alli = [alli; R.iti(:)]; %#ok<AGROW>
    end
end

wmean = @(M, w) (w(~isnan(w))' * M(~isnan(w), :)) ./ sum(w(~isnan(w)));
muAll = mean(psth, 1, "omitnan");
seAll = sem(psth);
mu13 = mean(psth(is13, :), 1, "omitnan");
mu18 = mean(psth(is18, :), 1, "omitnan");
muNw = wmean(psth, wN);
muTw = wmean(psth, wT);

% ----- fig1 pooled PSTH -----
figure("Position", [20, 20, 1600, 620], "Color", "w");
sgtitle({"Fig1  |  Pooled population PSTH  (all 18 movies)", ...
    "Black = mean of every movie  |  13 Aug t=0 is brightness peak  |  18 Aug t=0 is intensity dip"}, ...
    "FontSize", 13, "FontWeight", "bold");
subplot(1, 2, 1);
if any(is13)
    plot(tq, psth(is13, :)', "Color", [1.0 0.82 0.62], "LineWidth", 0.7); hold on;
end
if any(is18)
    plot(tq, psth(is18, :)', "Color", [0.70 0.88 0.90], "LineWidth", 0.7); hold on;
end
fill([tq fliplr(tq)], [muAll-seAll fliplr(muAll+seAll)], [0.2 0.2 0.2], ...
    "FaceAlpha", 0.18, "EdgeColor", "none");
p13 = plot(tq, mu13, "Color", col13, "LineWidth", 2.0);
p18 = plot(tq, mu18, "Color", col18, "LineWidth", 2.0);
pall = plot(tq, muAll, "k", "LineWidth", 2.8);
xline(0, "k-", "LineWidth", 1.4);
xlabel("Time from poke (s)"); ylabel("Population mean z");
title("Tint = movies  |  orange=0813  teal=0818  black=all 18");
box off;
try
    legend([p13 p18 pall], '0813 mean (12)', '0818 mean (6)', 'Pooled mean (18 movies)');
catch
end

subplot(1, 2, 2);
if any(is13)
    plot(tq, frac(is13, :)', "Color", [1.0 0.82 0.62], "LineWidth", 0.7); hold on;
end
if any(is18)
    plot(tq, frac(is18, :)', "Color", [0.70 0.88 0.90], "LineWidth", 0.7); hold on;
end
muF = mean(frac, 1, "omitnan"); seF = sem(frac);
fill([tq fliplr(tq)], [muF-seF fliplr(muF+seF)], [0.2 0.2 0.2], ...
    "FaceAlpha", 0.18, "EdgeColor", "none");
plot(tq, mean(frac(is13, :), 1, "omitnan"), "Color", col13, "LineWidth", 2.0);
plot(tq, mean(frac(is18, :), 1, "omitnan"), "Color", col18, "LineWidth", 2.0);
plot(tq, muF, "k", "LineWidth", 2.8);
xline(0, "k-", "LineWidth", 1.4);
xlabel("Time from poke (s)");
ylabel(sprintf("Fraction active (z>%g)", params.z_thresh));
title("Population recruitment");
box off;
exportgraphics(gcf, fullfile(avgDir, "fig1_population_psth.png"), "Resolution", 200);

% weighting inset
figure("Position", [20, 20, 980, 520], "Color", "w");
plot(tq, muAll, "k", "LineWidth", 2.6); hold on;
plot(tq, muNw, "Color", [0.45 0.45 0.45], "LineWidth", 1.8, "LineStyle", "--");
plot(tq, muTw, "Color", [0.20 0.20 0.20], "LineWidth", 1.6, "LineStyle", ":");
plot(tq, mu13, "Color", col13, "LineWidth", 1.8);
plot(tq, mu18, "Color", col18, "LineWidth", 1.8);
xline(0, "k-", "LineWidth", 1.3);
xlabel("Time from poke (s)"); ylabel("Population mean z");
title("Fig1b  |  Pooled mean is robust to weighting");
legend("Equal movie weight (primary)", "Neuron-weighted", "Trial-weighted", ...
    "0813 only", "0818 only", "Location", "best");
box off;
exportgraphics(gcf, fullfile(avgDir, "fig1b_weighting.png"), "Resolution", 200);

% ----- fig2 pooled STA -----
figure("Position", [20, 20, 980, 540], "Color", "w");
plot(tq, mean(staE, 1, "omitnan"), "Color", colExc, "LineWidth", 2.4); hold on;
plot(tq, mean(staI, 1, "omitnan"), "Color", colInh, "LineWidth", 2.4);
plot(tq, mean(staA, 1, "omitnan"), "k", "LineWidth", 2.0);
plot(tq, mean(staE(is13, :), 1, "omitnan"), "Color", colExc, "LineWidth", 1.0, "LineStyle", ":");
plot(tq, mean(staI(is13, :), 1, "omitnan"), "Color", colInh, "LineWidth", 1.0, "LineStyle", ":");
plot(tq, mean(staE(is18, :), 1, "omitnan"), "Color", colExc, "LineWidth", 1.0, "LineStyle", "--");
plot(tq, mean(staI(is18, :), 1, "omitnan"), "Color", colInh, "LineWidth", 1.0, "LineStyle", "--");
xline(0, "k-", "LineWidth", 1.4);
xlabel("Time from poke (s)"); ylabel("Mean STA (z)");
title(sprintf("Fig2  |  Pooled STA (mean of %d movies; dotted=0813, dashed=0818)", n));
legend("Excited (pooled)", "Inhibited (pooled)", "All cells (pooled)", ...
    "Location", "best");
box off;
exportgraphics(gcf, fullfile(avgDir, "fig2_neuron_sta.png"), "Resolution", 200);

% ----- fig3 pooled fractions -----
fePool = sum(nexc) / sum(nn);
fiPool = sum(ninh) / sum(nn);
figure("Position", [20, 20, 1500, 520], "Color", "w");
sgtitle(sprintf(['Fig3  |  Responsive cells pooled across %d movies  |  ', ...
    'neuron-pooled: %.1f%% exc, %.1f%% inh of %d cells'], ...
    n, 100*fePool, 100*fiPool, sum(nn)), "FontSize", 13, "FontWeight", "bold");
subplot(1, 3, 1);
b = bar([1 2], [mean(fe) mean(fi)]*100, 0.55, "EdgeColor", "k"); hold on;
b.FaceColor = "flat"; b.CData = [colExc; colInh];
errorbar([1 2], [mean(fe) mean(fi)]*100, [std(fe) std(fi)]*100, "k.", "LineWidth", 1.4);
plot([1 2], [100*fePool, 100*fiPool], "kd", "MarkerFaceColor", "y", "MarkerSize", 9);
set(gca, "XTick", 1:2, "XTickLabel", {"Excited", "Inhibited"});
ylabel("% of neurons");
title({"Bars = mean\pm SD of movies", "yellow diamond = all neurons pooled"});
box off;
subplot(1, 3, 2);
plot(100*fe(is13), 100*fi(is13), "o", "MarkerFaceColor", col13, ...
    "MarkerEdgeColor", "k", "MarkerSize", 8); hold on;
plot(100*fe(is18), 100*fi(is18), "s", "MarkerFaceColor", col18, ...
    "MarkerEdgeColor", "k", "MarkerSize", 8);
xlabel("% excited"); ylabel("% inhibited"); axis equal;
title("Each movie");
legend("0813", "0818", "Location", "best");
box off; grid on;
subplot(1, 3, 3);
histogram(log10(max(pp(is13), realmin)), "FaceColor", col13, "EdgeColor", "none"); hold on;
histogram(log10(max(pp(is18), realmin)), "FaceColor", col18, "EdgeColor", "none", "FaceAlpha", 0.7);
xlabel("log10 population p"); ylabel("# movies");
title(sprintf("%d / %d movies p<0.05", sum(pp<0.05), n));
box off;
exportgraphics(gcf, fullfile(avgDir, "fig3_responsive_cells.png"), "Resolution", 200);

% ----- fig4 aligned galleries (FOVs cannot be pixel-averaged) -----
g13 = "";
g18 = "";
if any(is13)
    g13 = fullfile(fileparts(alignFiles{find(is13, 1)}), "cohort_aligned_sig.png");
end
if any(is18)
    g18 = fullfile(fileparts(alignFiles{find(is18, 1)}), "cohort_aligned_sig.png");
end
figure("Position", [10, 10, 1700, 780], "Color", "w");
sgtitle({"Fig4  |  All recordings, aligned spatial maps", ...
    "Cannot average FOVs across mice; this is the pooled set. P left, A right."}, ...
    "FontSize", 13, "FontWeight", "bold");
subplot(1, 2, 1);
if exist(g13, "file") == 2
    image(imread(g13)); axis image off;
end
title(sprintf("13 Aug  |  %d movies  |  brightness peak", sum(is13)));
subplot(1, 2, 2);
if exist(g18, "file") == 2
    image(imread(g18)); axis image off;
end
title(sprintf("18 Aug  |  %d movies  |  intensity dip", sum(is18)));
exportgraphics(gcf, fullfile(avgDir, "fig4_spatial_aligned.png"), "Resolution", 160);
copyfile(fullfile(avgDir, "fig4_spatial_aligned.png"), fullfile(avgDir, "fig4_spatial_days.png"));

% ----- fig5 pooled timing -----
figure("Position", [20, 20, 1500, 520], "Color", "w");
sgtitle("Fig5  |  Pooled poke counts and intervals", "FontSize", 13, "FontWeight", "bold");
subplot(1, 2, 1);
bcol = repmat([0.6 0.6 0.6], n, 1);
bcol(is13, :) = repmat(col13, sum(is13), 1);
bcol(is18, :) = repmat(col18, sum(is18), 1);
bh = bar(1:n, npk, 0.82, "EdgeColor", "none");
bh.FaceColor = "flat";
bh.CData = bcol;
ylabel("# detected pokes"); xlabel("Recording (0813 then 0818)");
title(sprintf("Mean %.1f pokes / movie  |  total %d", mean(npk), sum(npk)));
box off;
subplot(1, 2, 2);
if ~isempty(alli)
    histogram(alli, "FaceColor", [0.35 0.40 0.70], "EdgeColor", "none");
    xlabel("Inter-poke interval (s)"); ylabel("Count");
    title(sprintf("Pooled ITI  median = %.2f s  (n=%d intervals)", median(alli), numel(alli)));
end
box off;
exportgraphics(gcf, fullfile(avgDir, "fig5_trial_timing.png"), "Resolution", 200);

fid = fopen(fullfile(avgDir, "pooled_numbers.txt"), "w");
fprintf(fid, "Pooled 13 Aug + 18 Aug peri-poke\n");
fprintf(fid, "movies=%d (0813=%d, 0818=%d)\n", n, sum(is13), sum(is18));
fprintf(fid, "neurons=%d\n", sum(nn));
fprintf(fid, "pokes=%d\n", sum(npk));
fprintf(fid, "movie-mean %%exc=%.2f  %%inh=%.2f\n", 100*mean(fe), 100*mean(fi));
fprintf(fid, "neuron-pooled %%exc=%.2f  %%inh=%.2f\n", 100*fePool, 100*fiPool);
fprintf(fid, "movies p<0.05: %d / %d\n", sum(pp<0.05), n);
fprintf(fid, "ITI median pooled=%.3f s\n", median(alli));
fprintf(fid, "pop z at t=0 (pooled mean)=%.4f\n", interp1(tq, muAll, 0));
fclose(fid);
fprintf("%s\n", fileread(fullfile(avgDir, "pooled_numbers.txt")));
close all;
end
