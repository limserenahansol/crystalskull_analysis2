function extract_stim_plot_aligned_maps(peakRoot, extRoot)
%EXTRACT_STIM_PLOT_ALIGNED_MAPS  Per-recording maps with mirror2 flipped.
%
% Aligns anterior to the RIGHT (mirror1 convention).

if nargin < 1 || isempty(peakRoot)
    peakRoot = fullfile(fileparts(mfilename("fullpath")), ...
        "stim_triggered_figures", "20250813_lab_meeting");
end
if nargin < 2 || isempty(extRoot)
    extRoot = "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\20250813";
end
idxFile = fullfile(peakRoot, "lab_meeting", "recording_index.csv");
outPer = fullfile(peakRoot, "lab_meeting", "aligned_per_recording");
outPair = fullfile(peakRoot, "lab_meeting", "aligned_pairs");
if exist(outPer, "dir") ~= 7
    mkdir(outPer);
end
if exist(outPair, "dir") ~= 7
    mkdir(outPair);
end

T = readtable(idxFile, "Delimiter", ",");
n = height(T);
cmapRb = red_blue_cmap();
sigMaps = cell(n, 1);
allMaps = cell(n, 1);
labels = cell(n, 1);

for i = 1:n
    movie = T.movie{i};
    label = T.label{i};
    labels{i} = label;
    flipLR = contains(label, "mirror2");
    sm = load_amp(peakRoot, movie);
    d = dir(fullfile(extRoot, sprintf("%s*_ext.mat", movie)));
    extFile = fullfile(extRoot, d(1).name);
    [Sn, frameSize] = extract_load_ext_spatial(extFile, numel(sm.amp), 2304);
    framePx = frameSize(1);
    sumimg = [];
    try
        sumimg = h5read(extFile, "/info/summary_image")';
    catch
    end
    respAll = reshape(full(Sn * sm.amp(:)), framePx, framePx)';
    ampSig = sm.amp(:);
    ampSig(~(sm.isExc(:) | sm.isInh(:))) = 0;
    respSig = reshape(full(Sn * ampSig), framePx, framePx)';
    if flipLR
        respAll = fliplr(respAll);
        respSig = fliplr(respSig);
        if ~isempty(sumimg)
            sumimg = fliplr(sumimg);
        end
    end
    allMaps{i} = imresize(respAll, [256, 256]);
    sigMaps{i} = imresize(respSig, [256, 256]);
    write_three_panel(sumimg, respAll, respSig, sm, label, flipLR, cmapRb, ...
        fullfile(outPer, sprintf("%s_fig4_aligned.png", label)));
    fprintf("  %s  flipLR=%d  exc=%d inh=%d\n", label, flipLR, sum(sm.isExc), sum(sm.isInh));
end

write_cohort(sigMaps, labels, cmapRb, fullfile(outPer, "cohort_aligned_sig.png"));
write_pairs(T, sigMaps, cmapRb, outPair);
fprintf("Aligned maps -> %s\n", outPer);
end

function write_three_panel(sumimg, respAll, respSig, sm, label, flipLR, cmapRb, outFile)
showPx = 768;
respAllShow = imresize(respAll, [showPx, showPx]);
respSigShow = imresize(respSig, [showPx, showPx]);
sigAbs = abs(respSigShow(respSigShow ~= 0));
if isempty(sigAbs)
    sigAbs = abs(respAllShow(respAllShow ~= 0));
end
climv = prctile(sigAbs, 95);
if isempty(climv) || ~isfinite(climv) || climv == 0
    climv = 1;
end
    alignNote = "aligned: posterior left, anterior right";
    if flipLR
        alignNote = "aligned: posterior left, anterior right  (mirror2 flipped)";
    end
figure("Position", [20, 20, 1680, 580], "Color", "w");
sgtitle({label, sprintf("Spatial poke map  |  %s  |  red=exc  blue=inh", alignNote)}, ...
    "FontSize", 12, "FontWeight", "bold", "Interpreter", "none");
subplot(1, 3, 1);
if ~isempty(sumimg)
    imagesc(imresize(double(sumimg), [showPx, showPx]));
    colormap(gca, "gray");
    axis image off;
else
    axis off;
end
title("Anatomy");
add_ap_labels();
subplot(1, 3, 2);
imagesc(respAllShow);
colormap(gca, cmapRb);
clim([-climv, climv]);
axis image off;
cb = colorbar;
cb.Label.String = "amp (post-pre, z)";
title(sprintf("All neurons (n=%d)", numel(sm.amp)));
add_ap_labels();
subplot(1, 3, 3);
imagesc(respSigShow);
colormap(gca, cmapRb);
clim([-climv, climv]);
axis image off;
colorbar;
title(sprintf("Significant  exc=%d  inh=%d", sum(sm.isExc), sum(sm.isInh)));
add_ap_labels();
exportgraphics(gcf, outFile, "Resolution", 200);
close(gcf);
end

function add_ap_labels()
text(0.04, 0.96, "P", "Units", "normalized", "Color", [0.1 0.35 0.85], ...
    "FontWeight", "bold", "FontSize", 12, "VerticalAlignment", "top");
text(0.96, 0.96, "A", "Units", "normalized", "Color", [0.85 0.20 0.15], ...
    "FontWeight", "bold", "FontSize", 12, "VerticalAlignment", "top", ...
    "HorizontalAlignment", "right");
end

function write_cohort(sigMaps, labels, cmapRb, outFile)
n = numel(sigMaps);
ncol = 4;
nrow = ceil(n / ncol);
figure("Position", [20, 20, 1680, 980], "Color", "w");
sgtitle({"All recordings  |  aligned poke maps (significant cells)", ...
    "mirror2 flipped so anterior is on the RIGHT in every panel   |   red=excited  blue=inhibited"}, ...
    "FontSize", 13, "FontWeight", "bold");
for i = 1:n
    subplot(nrow, ncol, i);
    M = sigMaps{i};
    climv = prctile(abs(M(M ~= 0)), 95);
    if isempty(climv) || ~isfinite(climv) || climv == 0
        climv = 1;
    end
    imagesc(M);
    colormap(gca, cmapRb);
    clim([-climv, climv]);
    axis image off;
    add_ap_labels();
    title(strrep(labels{i}, "_", "\_"), "FontSize", 8);
end
exportgraphics(gcf, outFile, "Resolution", 180);
close(gcf);
end

function write_pairs(T, sigMaps, cmapRb, outPair)
keys = unique(regexprep(T.label, "_mirror[12]$", ""));
for k = 1:numel(keys)
    key = keys{k};
    i1 = find(strcmp(T.label, sprintf("%s_mirror1", key)), 1);
    i2 = find(strcmp(T.label, sprintf("%s_mirror2", key)), 1);
    if isempty(i1) || isempty(i2)
        continue
    end
    figure("Position", [20, 20, 1200, 520], "Color", "w");
    sgtitle({strrep(key, "_", " "), ...
        "Same session, two optical paths  |  mirror2 flipped to match  |  P left, A right"}, ...
        "FontSize", 12, "FontWeight", "bold");
    subplot(1, 2, 1);
    show_sig(sigMaps{i1}, cmapRb, "mirror1 (native)");
    subplot(1, 2, 2);
    show_sig(sigMaps{i2}, cmapRb, "mirror2 (flipped)");
    exportgraphics(gcf, fullfile(outPair, sprintf("%s_mirrors.png", key)), "Resolution", 200);
    close(gcf);
end
end

function show_sig(M, cmapRb, ttl)
climv = prctile(abs(M(M ~= 0)), 95);
if isempty(climv) || ~isfinite(climv) || climv == 0
    climv = 1;
end
imagesc(M);
colormap(gca, cmapRb);
clim([-climv, climv]);
axis image off;
title(ttl);
add_ap_labels();
end

function cmapRb = red_blue_cmap()
halfN = 128;
cmapRb = [[linspace(0, 1, halfN), ones(1, halfN)]' ...
    [linspace(0, 1, halfN), linspace(1, 0, halfN)]' ...
    [ones(1, halfN), linspace(1, 0, halfN)]'];
end

function sm = load_amp(peakRoot, movie)
pMap = fullfile(peakRoot, movie, "spatial_map.mat");
pRes = fullfile(peakRoot, movie, "stim_results.mat");
if exist(pMap, "file") == 2
    sm = load(pMap);
elseif exist(pRes, "file") == 2
    tmp = load(pRes, "R");
    sm.amp = tmp.R.amp;
    sm.pval = tmp.R.pval;
else
    error("extract_stim_plot_aligned_maps:NoAmp", ...
        "Need spatial_map.mat or stim_results.mat in %s", movie);
end
if ~isfield(sm, "isExc") || ~isfield(sm, "isInh")
    sm.isExc = sm.pval < 0.05 & sm.amp > 0;
    sm.isInh = sm.pval < 0.05 & sm.amp < 0;
end
end
