function extract_stim_ap_takehome()
%EXTRACT_STIM_AP_TAKEHOME  One lab figure: anterior vs posterior after alignment.

addpath(fileparts(mfilename("fullpath")));
root = fileparts(mfilename("fullpath"));
outDir = fullfile(root, "stim_triggered_figures", "combined_0813_0818", "n18");
if exist(outDir, "dir") ~= 7
    mkdir(outDir);
end

jobs = {
    fullfile(root, "stim_triggered_figures", "20250813_lab_meeting"), ...
        "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\20250813"
    fullfile(root, "stim_triggered_figures", "20250818_lab_meeting"), ...
        "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\20250818"
    };

maps = {};
ampP = [];
ampA = [];
n = 0;
for j = 1:size(jobs, 1)
    peakRoot = jobs{j, 1};
    extRoot = jobs{j, 2};
    T = readtable(fullfile(peakRoot, "lab_meeting", "recording_index.csv"), "Delimiter", ",");
    for i = 1:height(T)
        movie = char(T.movie{i});
        label = char(T.label{i});
        tmp = load(fullfile(peakRoot, movie, "stim_results.mat"), "R");
        amp = tmp.R.amp(:);
        pval = tmp.R.pval(:);
        isResp = pval < 0.05;
        d = dir(fullfile(extRoot, sprintf("%s*_ext.mat", movie)));
        if isempty(d)
            warning("No EXTRACT file for %s", movie);
            continue
        end
        [Sn, frameSize] = extract_load_ext_spatial(fullfile(extRoot, d(1).name), numel(amp), 2304);
        H = frameSize(1);
        nPix = size(Sn, 1);
        rowS = rem((0:nPix-1)', H) + 1;
        flipLR = contains(label, "mirror2");
        colSum = full(sum(Sn, 1));
        if flipLR
            cx = full(((H + 1) - rowS)' * Sn) ./ max(colSum, eps);
        else
            cx = full((rowS' * Sn) ./ max(colSum, eps));
        end
        ok = colSum(:) > 0 & isfinite(amp) & isfinite(cx(:));
        cx = cx(:);
        q = quantile(cx(ok), [1/3, 2/3]);
        isP = ok & cx <= q(1);
        isA = ok & cx >= q(2);
        n = n + 1;
        ampP(n, 1) = mean(amp(isP), "omitnan"); %#ok<AGROW>
        ampA(n, 1) = mean(amp(isA), "omitnan"); %#ok<AGROW>
        ampSig = amp;
        ampSig(~isResp) = 0;
        M = reshape(full(Sn * ampSig), H, H)';
        if flipLR
            M = fliplr(M);
        end
        maps{n} = imresize(M, [192, 192]); %#ok<AGROW>
        fprintf("  %s  P=%.4f  A=%.4f  d=%.4f  nP=%d nA=%d\n", ...
            label, ampP(n), ampA(n), ampA(n)-ampP(n), sum(isP), sum(isA));
    end
end

dlt = ampA - ampP;
p = signrank(dlt);
nUp = sum(dlt > 0);
fprintf("n=%d  A>P in %d/%d  median(A-P)=%.4f  signrank p=%.3g\n", ...
    n, nUp, n, median(dlt), p);

stack = zeros(size(maps{1}));
for i = 1:n
    M = maps{i};
    s = prctile(abs(M(M ~= 0)), 95);
    if ~isfinite(s) || s == 0
        s = 1;
    end
    stack = stack + M ./ s;
end
meanMap = stack / n;

colExc = [0.85 0.20 0.15];
colInh = [0.10 0.35 0.85];
halfN = 128;
cmapRb = [[linspace(0, 1, halfN), ones(1, halfN)]' ...
    [linspace(0, 1, halfN), linspace(1, 0, halfN)]' ...
    [ones(1, halfN), linspace(1, 0, halfN)]'];

figure("Position", [40, 40, 1680, 720], "Color", "w");
sgtitle({"Take-home:  anterior is more excited,  posterior is more inhibited", ...
    sprintf("%d / %d recordings  anterior > posterior   |   signrank p = %.2g", nUp, n, p)}, ...
    "FontSize", 18, "FontWeight", "bold");

subplot(1, 2, 1);
imagesc(meanMap);
colormap(gca, cmapRb);
c = prctile(abs(meanMap(:)), 99);
if ~isfinite(c) || c == 0; c = 1; end
clim([-c, c]);
axis image off;
title("Mean map after alignment   (P left, A right)", "FontSize", 16);
text(0.04, 0.96, "P", "Units", "normalized", "Color", colInh, ...
    "FontWeight", "bold", "FontSize", 28, "VerticalAlignment", "top");
text(0.96, 0.96, "A", "Units", "normalized", "Color", colExc, ...
    "FontWeight", "bold", "FontSize", 28, "VerticalAlignment", "top", ...
    "HorizontalAlignment", "right");
cb = colorbar;
cb.Label.String = "poke amplitude";
cb.FontSize = 12;

subplot(1, 2, 2);
plot([1 2], [ampP ampA]', "-", "Color", [0.72 0.72 0.72], "LineWidth", 1.4); hold on;
plot(ones(n,1), ampP, "o", "MarkerFaceColor", colInh, "MarkerEdgeColor", "k", "MarkerSize", 9);
plot(2*ones(n,1), ampA, "o", "MarkerFaceColor", colExc, "MarkerEdgeColor", "k", "MarkerSize", 9);
plot([1 2], [mean(ampP) mean(ampA)], "k-", "LineWidth", 3.2);
plot(1, mean(ampP), "kd", "MarkerFaceColor", "k", "MarkerSize", 12);
plot(2, mean(ampA), "kd", "MarkerFaceColor", "k", "MarkerSize", 12);
yline(0, "k:", "LineWidth", 1);
xlim([0.6 2.4]);
set(gca, "XTick", [1 2], "XTickLabel", {"Posterior", "Anterior"}, "FontSize", 16);
ylabel("Mean poke amplitude (z)", "FontSize", 16);
title({"Each gray line = one recording", "Cells in the back third vs front third of the FOV"}, "FontSize", 14);
box off;

outFile = fullfile(outDir, "takehome_anterior_posterior.png");
exportgraphics(gcf, outFile, "Resolution", 220);
try
    copyfile(outFile, fullfile(root, "stim_triggered_figures", "20250813_lab_meeting", ...
        "lab_meeting", "cohort_average", "takehome_anterior_posterior.png"));
catch
end
close(gcf);
fprintf("Wrote %s\n", outFile);
end
