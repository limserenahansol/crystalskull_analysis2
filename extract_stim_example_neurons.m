function extract_stim_example_neurons()
%EXTRACT_STIM_EXAMPLE_NEURONS  Single-cell traces + FOV maps for Mark.
%
% For each 13 Aug and 18 Aug recording, pick two excited and two inhibited
% neurons (strongest in the anterior half and posterior half), then write:
%   per_recording/<label>/example_traces.png
%   per_recording/<label>/example_maps.png
%   per_recording/<label>/example_fov.png
%   all_example_neurons.csv

root = fileparts(mfilename("fullpath"));
addpath(root);

outDir = fullfile(root, "stim_triggered_figures", "combined_0813_0818", "example_neurons");
perDir = fullfile(outDir, "per_recording");
if exist(perDir, "dir") ~= 7
    mkdir(perDir);
end

jobs = {
    fullfile(root, "stim_triggered_figures", "20250813_lab_meeting"), ...
        "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\20250813", ...
        "13 Aug"
    fullfile(root, "stim_triggered_figures", "20250818_lab_meeting"), ...
        "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\20250818", ...
        "18 Aug"
    };

masterFile = fullfile(outDir, "all_example_neurons.csv");
fid = fopen(masterFile, "w");
fprintf(fid, "day,label,movie,tag,neuron_id,kind,amp,centroid_x,n_exc,n_inh,p_pop\n");

nOk = 0;
nAll = 0;
allLabels = {};
allDays = {};
for j = 1:size(jobs, 1)
    peakRoot = jobs{j, 1};
    extRoot = jobs{j, 2};
    dayLabel = jobs{j, 3};
    if strcmp(dayLabel, "13 Aug")
        stimMat = fullfile(extRoot, "estimated_stim_timing_20250813.mat");
    else
        stimMat = fullfile(extRoot, "estimated_stim_timing_20250818.mat");
    end
    T = readtable(fullfile(peakRoot, "lab_meeting", "recording_index.csv"), "Delimiter", ",");
    T.dayLabel = repmat({dayLabel}, height(T), 1);
    St = load_stim_file(stimMat);
    nAll = nAll + height(T);

    for i = 1:height(T)
        movie = char(T.movie{i});
        recLabel = char(T.label{i});
        recOut = fullfile(perDir, recLabel);
        if exist(recOut, "dir") ~= 7
            mkdir(recOut);
        end
        dispLabel = recording_title(T, i);
        fprintf("\n==== %s  %d/%d  %s ====\n", dayLabel, i, height(T), recLabel);
        try
            tracesFile = fullfile(recOut, "example_traces.png");
            csvFile = fullfile(recOut, "example_neurons.csv");
            if exist(tracesFile, "file") == 2 && exist(csvFile, "file") == 2
                fprintf("  traces exist; remaking high-contrast FOV maps\n");
                examples = examples_from_csv(csvFile, peakRoot, movie);
                remap_maps(peakRoot, extRoot, movie, recLabel, dispLabel, recOut, examples);
            else
                examples = process_one(peakRoot, extRoot, St, movie, recLabel, dispLabel, recOut);
            end
            nExc = examples(1).nExc;
            nInh = examples(1).nInh;
            pPop = examples(1).pPop;
            for k = 1:numel(examples)
                e = examples(k);
                fprintf(fid, "%s,%s,%s,%s,%d,%s,%.5f,%.2f,%d,%d,%.3e\n", ...
                    dayLabel, recLabel, movie, char(e.tag), e.id, e.kind, e.amp, e.cx, nExc, nInh, pPop);
            end
            allLabels{end+1} = recLabel; %#ok<AGROW>
            allDays{end+1} = dayLabel; %#ok<AGROW>
            nOk = nOk + 1;
        catch ME
            warning("extract_stim_example_neurons:MovieFailed", "%s: %s", recLabel, ME.message);
        end
    end
end
fclose(fid);

write_overview(perDir, allLabels, allDays, outDir);
fprintf("\nDone  %d/%d recordings  -> %s\n", nOk, nAll, outDir);
end

function St = load_stim_file(stimMat)
raw = load(stimMat);
St.fnames = cellfun(@char, raw.fnames, "UniformOutput", false);
St.peakIdx = cellfun(@(v) double(v(:)), raw.peak_idx, "UniformOutput", false);
end

function txt = recording_title(T, i)
mouse = char(string(T.mouse{i}));
session = char(string(T.session{i}));
mirror = char(string(T.mirror{i}));
pwrTxt = power_label(T, i);
dayLabel = "13 Aug";
if ismember("dayLabel", T.Properties.VariableNames)
    dayLabel = char(string(T.dayLabel{i}));
end
txt = sprintf("%s  |  %s  |  %s  |  %s  |  %s", mouse, dayLabel, session, pwrTxt, mirror);
end

function pwrTxt = power_label(T, i)
raw = T.power_mWcm2;
if iscell(raw)
    val = raw{i};
else
    val = raw(i);
end
sVal = string(val);
if any(ismissing(sVal)) || strlength(strtrim(sVal(1))) == 0
    pwrTxt = "pwr NA";
    return
end
s = strtrim(char(sVal(1)));
if any(strcmpi(s, ["unknown", "nan", "NA", "<missing>"]))
    pwrTxt = "pwr NA";
    return
end
num = str2double(s);
if isfinite(num)
    pwrTxt = sprintf("%g mW/cm2", num);
else
    pwrTxt = sprintf("%s mW/cm2", s);
end
end

function examples = process_one(peakRoot, extRoot, St, movie, recLabel, dispLabel, recOut)
res = load(fullfile(peakRoot, movie, "stim_results.mat"), "R");
R = res.R;
amp = R.amp(:);
pval = R.pval(:);
isExc = pval < 0.05 & amp > 0;
isInh = pval < 0.05 & amp < 0;
flipLR = contains(recLabel, "mirror2");

extFile = char(R.ext_file);
if exist(extFile, "file") ~= 2
    d = dir(fullfile(extRoot, sprintf("%s*_ext.mat", movie)));
    if isempty(d)
        error("extract_stim_example_neurons:NoExt", "EXTRACT file not found for %s.", movie);
    end
    extFile = fullfile(extRoot, d(1).name);
end

onset = poke_frames_for_movie(St, movie);
fprintf("  Loading EXTRACT traces...\n");
tw = double(h5read(extFile, "/temporal_weights"));
if size(tw, 1) < size(tw, 2) && size(tw, 2) >= max(onset)
    tw = tw.';
end
nFrames = size(tw, 1);
nNeurons = size(tw, 2);
if nNeurons ~= numel(amp)
    error("extract_stim_example_neurons:CountMismatch", ...
        "%s: traces have %d neurons, results have %d.", movie, nNeurons, numel(amp));
end
fs = 33;

sumimg = [];
try
    sumimg = h5read(extFile, "/info/summary_image")';
catch
end
if flipLR && ~isempty(sumimg)
    sumimg = fliplr(sumimg);
end

fprintf("  Reconstructing footprints...\n");
[Sn, frameSize] = extract_load_ext_spatial(extFile, nNeurons, 2304);
H = frameSize(1);
nPix = size(Sn, 1);
rowS = rem((0:nPix-1)', H) + 1;
colSum = full(sum(Sn, 1));
if flipLR
    cx = full(((H + 1) - rowS)' * Sn) ./ max(colSum, eps);
else
    cx = full((rowS' * Sn) ./ max(colSum, eps));
end
cx = cx(:);
ok = colSum(:) > 0 & isfinite(amp) & isfinite(cx);

examples = pick_examples(amp, isExc, isInh, cx, ok);
for k = 1:numel(examples)
    examples(k).nExc = sum(isExc);
    examples(k).nInh = sum(isInh);
    examples(k).pPop = R.p_pop;
end
write_example_table(examples, recOut);

ids = [examples.id];
zSel = (tw(:, ids) - mean(tw(:, ids), 1)) ./ max(std(tw(:, ids), 0, 1), 1e-8);
periWin = [-2, 3];
pk = round(periWin(1)*fs):round(periWin(2)*fs);
tPeri = pk(:) / fs;
sta = compute_sta_cells(zSel, onset, pk, nFrames);

write_trace_figure(tw, onset, fs, sta, tPeri, examples, dispLabel, recOut);
write_map_figure(sumimg, Sn, H, amp, isExc, isInh, examples, dispLabel, recOut, flipLR);
end

function remap_maps(peakRoot, extRoot, movie, recLabel, dispLabel, recOut, examples)
res = load(fullfile(peakRoot, movie, "stim_results.mat"), "R");
R = res.R;
amp = R.amp(:);
pval = R.pval(:);
isExc = pval < 0.05 & amp > 0;
isInh = pval < 0.05 & amp < 0;
flipLR = contains(recLabel, "mirror2");

extFile = char(R.ext_file);
if exist(extFile, "file") ~= 2
    d = dir(fullfile(extRoot, sprintf("%s*_ext.mat", movie)));
    if isempty(d)
        error("extract_stim_example_neurons:NoExt", "EXTRACT file not found for %s.", movie);
    end
    extFile = fullfile(extRoot, d(1).name);
end

sumimg = [];
try
    sumimg = h5read(extFile, "/info/summary_image")';
catch
end
if flipLR && ~isempty(sumimg)
    sumimg = fliplr(sumimg);
end

fprintf("  Reconstructing footprints for map...\n");
[Sn, frameSize] = extract_load_ext_spatial(extFile, numel(amp), 2304);
write_map_figure(sumimg, Sn, frameSize(1), amp, isExc, isInh, examples, dispLabel, recOut, flipLR);
end

function onset = poke_frames_for_movie(St, movie)
movie = char(movie);
hit = 0;
for i = 1:numel(St.fnames)
    stem = regexprep(St.fnames{i}, "\.h5$|\.mat$", "");
    if strcmp(stem, movie) || startsWith(stem, movie) || startsWith(movie, stem)
        hit = i;
        break
    end
end
if hit == 0
    error("extract_stim_example_neurons:NoStim", "No stim entry for %s.", movie);
end
onset = round(St.peakIdx{hit});
onset = onset(isfinite(onset) & onset >= 1);
end

function examples = examples_from_csv(csvFile, peakRoot, movie)
T = readtable(csvFile, "Delimiter", ",");
res = load(fullfile(peakRoot, movie, "stim_results.mat"), "R");
R = res.R;
colExc = [0.85, 0.20, 0.15];
colInh = [0.10, 0.35, 0.85];
examples = struct("id", {}, "kind", {}, "tag", {}, "amp", {}, "cx", {}, "color", {}, ...
    "nExc", {}, "nInh", {}, "pPop", {});
nExc = sum(R.pval(:) < 0.05 & R.amp(:) > 0);
nInh = sum(R.pval(:) < 0.05 & R.amp(:) < 0);
for i = 1:height(T)
    kind = char(string(T.kind{i}));
    if strcmp(kind, "excited")
        color = colExc;
    else
        color = colInh;
    end
    examples(end+1) = pack_example(T.neuron_id(i), kind, string(T.tag{i}), ...
        T.amp(i), T.centroid_x(i), color); %#ok<AGROW>
    examples(end).nExc = nExc;
    examples(end).nInh = nInh;
    examples(end).pPop = R.p_pop;
end
end

function examples = pick_examples(amp, isExc, isInh, cx, ok)
examples = struct("id", {}, "kind", {}, "tag", {}, "amp", {}, "cx", {}, "color", {}, ...
    "nExc", {}, "nInh", {}, "pPop", {});
colExc = [0.85, 0.20, 0.15];
colInh = [0.10, 0.35, 0.85];
examples = add_pair(examples, amp, isExc & ok, cx, "excited", colExc, true);
examples = add_pair(examples, amp, isInh & ok, cx, "inhibited", colInh, false);
end

function examples = add_pair(examples, amp, mask, cx, kind, color, wantPositive)
ids = find(mask);
if numel(ids) < 2
    error("extract_stim_example_neurons:TooFew", "Need at least 2 %s cells.", kind);
end
medX = median(cx(ids));
if wantPositive
    score = amp;
else
    score = -amp;
end
antIds = ids(cx(ids) >= medX);
postIds = ids(cx(ids) < medX);
if isempty(antIds)
    antIds = ids;
end
if isempty(postIds)
    postIds = ids;
end
[~, ia] = max(score(antIds));
[~, ip] = max(score(postIds));
idAnt = antIds(ia);
idPost = postIds(ip);
if idAnt == idPost
    [~, ord] = sort(score(ids), "descend");
    idAnt = ids(ord(1));
    idPost = ids(ord(min(2, numel(ord))));
end
if strcmp(kind, "excited")
    tagAnt = "E1 anterior";
    tagPost = "E2 posterior";
else
    tagAnt = "I1 anterior";
    tagPost = "I2 posterior";
end
examples(end+1) = pack_example(idAnt, kind, tagAnt, amp(idAnt), cx(idAnt), color);
examples(end+1) = pack_example(idPost, kind, tagPost, amp(idPost), cx(idPost), color);
end

function s = pack_example(id, kind, tag, amp, cx, color)
s = struct("id", id, "kind", kind, "tag", tag, "amp", amp, "cx", cx, "color", color, ...
    "nExc", 0, "nInh", 0, "pPop", NaN);
end

function sta = compute_sta_cells(zSel, onset, pk, nFrames)
nEx = size(zSel, 2);
nP = numel(pk);
sta = zeros(nEx, nP);
nT = 0;
for k = 1:numel(onset)
    idx = onset(k) + pk;
    if idx(1) < 1 || idx(end) > nFrames
        continue
    end
    sta = sta + zSel(idx, :).';
    nT = nT + 1;
end
if nT > 0
    sta = sta / nT;
end
end

function write_example_table(examples, outDir)
fid = fopen(fullfile(outDir, "example_neurons.csv"), "w");
fprintf(fid, "tag,neuron_id,kind,amp,centroid_x\n");
for i = 1:numel(examples)
    e = examples(i);
    fprintf(fid, "%s,%d,%s,%.5f,%.2f\n", char(e.tag), e.id, e.kind, e.amp, e.cx);
end
fclose(fid);
end

function write_trace_figure(tw, onset, fs, sta, tPeri, examples, label, outDir)
nEx = numel(examples);
t = (0:size(tw, 1)-1)' / fs;
winSec = 45;
tStart = choose_trace_window(onset, fs, winSec, size(tw, 1));
tEnd = tStart + winSec;
inWin = t >= tStart & t < tEnd;
pokeWin = onset(onset / fs >= tStart & onset / fs < tEnd) / fs;

figure("Position", [20, 20, 1760, 980], "Color", "w", "Visible", "off");
sgtitle({sprintf("Example calcium traces   |   %s", label), ...
    "EXTRACT activity with poke times (left) and peri-poke STA of the same cell (right)"}, ...
    "FontSize", 14, "FontWeight", "bold", "Interpreter", "none");

for i = 1:nEx
    e = examples(i);
    tr = tw(inWin, e.id);
    tr = tr - median(tr, "omitnan");
    ax = subplot(nEx, 2, 2*i-1);
    hold(ax, "on");
    plot(ax, t(inWin), tr, "Color", e.color, "LineWidth", 1.15);
    if ~isempty(pokeWin)
        xline(ax, pokeWin, "-", "Color", [0.55, 0.55, 0.55], "LineWidth", 0.7);
    end
    ylabel(ax, "dF/F-like");
    title(ax, sprintf("%s   |   neuron %d   |   amp = %+.3f", char(e.tag), e.id, e.amp), ...
        "Interpreter", "none", "Color", e.color, "FontSize", 11);
    box(ax, "off");
    xlim(ax, [tStart, tEnd]);
    if i == nEx
        xlabel(ax, "Time (s)");
    else
        set(ax, "XTickLabel", []);
    end
    if i == 1
        text(ax, 0.02, 0.92, "vertical lines = pokes", "Units", "normalized", ...
            "Color", [0.35, 0.35, 0.35], "FontSize", 9);
    end

    ax2 = subplot(nEx, 2, 2*i);
    hold(ax2, "on");
    yline(ax2, 0, "k:", "LineWidth", 0.8);
    xline(ax2, 0, "k-", "LineWidth", 1.2);
    plot(ax2, tPeri, sta(i, :), "Color", e.color, "LineWidth", 2.2);
    ylabel(ax2, "z");
    title(ax2, "Peri-poke STA", "FontSize", 11);
    box(ax2, "off");
    xlim(ax2, [tPeri(1), tPeri(end)]);
    if i == nEx
        xlabel(ax2, "Time from poke (s)");
    else
        set(ax2, "XTickLabel", []);
    end
end

exportgraphics(gcf, fullfile(outDir, "example_traces.png"), "Resolution", 200);
close(gcf);
end

function tStart = choose_trace_window(onset, fs, winSec, nFrames)
tPoke = onset / fs;
tMax = nFrames / fs;
best = 1;
bestN = -1;
for i = 1:numel(tPoke)
    t0 = max(tPoke(i) - 2, 0);
    if t0 + winSec > tMax
        break
    end
    nHere = sum(tPoke >= t0 & tPoke < t0 + winSec);
    if nHere > bestN
        bestN = nHere;
        best = t0;
    end
end
tStart = best;
end

function write_map_figure(sumimg, Sn, H, amp, isExc, isInh, examples, label, outDir, flipLR)
showPx = 768;
halfN = 128;
cmapRb = [[linspace(0, 1, halfN), ones(1, halfN)]' ...
    [linspace(0, 1, halfN), linspace(1, 0, halfN)]' ...
    [ones(1, halfN), linspace(1, 0, halfN)]'];

ampSig = amp;
ampSig(~(isExc | isInh)) = 0;
respSig = reshape(full(Sn * ampSig), H, H)';
if flipLR
    respSig = fliplr(respSig);
end
respShow = imresize(respSig, [showPx, showPx]);
sigAbs = abs(respShow(respShow ~= 0));
climv = prctile(sigAbs, 95);
if isempty(climv) || ~isfinite(climv) || climv == 0
    climv = 1;
end

figure("Position", [20, 20, 1760, 820], "Color", "w", "Visible", "off");
sgtitle({sprintf("Where these neurons reside   |   %s", label), ...
    "Crystal-skull FOV (mesoscope summary image). P = posterior, A = anterior."}, ...
    "FontSize", 14, "FontWeight", "bold", "Interpreter", "none");

ax1 = subplot(1, 2, 1);
show_fov(ax1, sumimg, showPx);
hold(ax1, "on");
overlay_examples(ax1, Sn, H, showPx, examples, true, flipLR);
title(ax1, "Crystal-skull FOV  +  example neurons", "FontSize", 13);
add_ap_labels(ax1);

ax2 = subplot(1, 2, 2);
imagesc(respShow);
colormap(ax2, cmapRb);
clim(ax2, [-climv, climv]);
axis(ax2, "image");
axis(ax2, "off");
hold(ax2, "on");
overlay_examples(ax2, Sn, H, showPx, examples, false, flipLR);
title(ax2, sprintf("All significant cells   exc=%d   inh=%d", sum(isExc), sum(isInh)), ...
    "FontSize", 13);
add_ap_labels(ax2);
cb = colorbar(ax2);
cb.Label.String = "amp (post - pre, z)";

exportgraphics(gcf, fullfile(outDir, "example_maps.png"), "Resolution", 200);
close(gcf);

figFov = figure("Position", [20, 20, 820, 780], "Color", "w", "Visible", "off");
axF = axes(figFov);
show_fov(axF, sumimg, showPx);
hold(axF, "on");
overlay_examples(axF, Sn, H, showPx, examples, true, flipLR);
title(axF, label, "FontSize", 11, "Interpreter", "none");
add_ap_labels(axF);
exportgraphics(figFov, fullfile(outDir, "example_fov.png"), "Resolution", 180);
close(figFov);
end

function show_fov(ax, sumimg, showPx)
if isempty(sumimg)
    imagesc(ax, zeros(showPx, showPx));
else
    imagesc(ax, enhance_fov(sumimg, showPx));
end
colormap(ax, "gray");
clim(ax, [0, 1]);
axis(ax, "image");
axis(ax, "off");
end

function img = enhance_fov(sumimg, showPx)
img = imresize(double(sumimg), [showPx, showPx]);
vals = img(isfinite(img));
if isempty(vals)
    img = zeros(showPx, showPx);
    return
end
lo = prctile(vals, 2);
hi = prctile(vals, 97);
if ~(isfinite(lo) && isfinite(hi)) || hi <= lo
    lo = min(vals);
    hi = max(vals);
end
img = (img - lo) / max(hi - lo, eps);
img = min(max(img, 0), 1);
img = img .^ 0.5;
img = adapthisteq(img, "NumTiles", [16, 16], "ClipLimit", 0.04, "Range", "original");
img = min(max(img, 0), 1);
end

function overlay_examples(ax, Sn, H, showPx, examples, fillFoot, flipLR)
scale = showPx / H;
for i = 1:numel(examples)
    e = examples(i);
    w = reshape(full(Sn(:, e.id)), H, H)';
    if flipLR
        w = fliplr(w);
    end
    w = imresize(w, [showPx, showPx]);
    [rr, cc] = find(w > 0);
    if isempty(rr)
        continue
    end
    ww = w(w > 0);
    cy = sum(rr .* ww) / sum(ww);
    cx = sum(cc .* ww) / sum(ww);
    if fillFoot
        mask = w >= 0.15 * max(w(:));
        [B, ~] = bwboundaries(mask, "noholes");
        for b = 1:numel(B)
            plot(ax, B{b}(:, 2), B{b}(:, 1), "Color", e.color, "LineWidth", 1.8);
        end
    end
    plot(ax, cx, cy, "o", "MarkerSize", 11, "LineWidth", 1.6, ...
        "MarkerEdgeColor", "k", "MarkerFaceColor", e.color);
    text(ax, cx + 14*scale, cy, char(e.tag), "Color", e.color, "FontWeight", "bold", ...
        "FontSize", 11, "Interpreter", "none", "BackgroundColor", [1, 1, 1]);
end
end

function add_ap_labels(ax)
text(ax, 0.04, 0.96, "P", "Units", "normalized", "Color", [0.10, 0.35, 0.85], ...
    "FontWeight", "bold", "FontSize", 20, "VerticalAlignment", "top");
text(ax, 0.96, 0.96, "A", "Units", "normalized", "Color", [0.85, 0.20, 0.15], ...
    "FontWeight", "bold", "FontSize", 20, "VerticalAlignment", "top", ...
    "HorizontalAlignment", "right");
end

function write_overview(perDir, allLabels, allDays, outDir)
n = numel(allLabels);
ncol = 6;
nrow = ceil(n / ncol);
figure("Position", [20, 20, 2000, 1100], "Color", "w", "Visible", "off");
sgtitle("Example neurons on each crystal-skull FOV   |   13 Aug + 18 Aug, 18 recordings", ...
    "FontSize", 16, "FontWeight", "bold");
for i = 1:n
    recLabel = allLabels{i};
    fovFile = fullfile(perDir, recLabel, "example_fov.png");
    ax = subplot(nrow, ncol, i);
    if exist(fovFile, "file") == 2
        imshow(imread(fovFile));
    else
        axis(ax, "off");
    end
    title(ax, sprintf("%s  %s", allDays{i}, strrep(recLabel, "_", " ")), ...
        "FontSize", 8, "Interpreter", "none");
end
exportgraphics(gcf, fullfile(outDir, "example_fov_overview_n18.png"), "Resolution", 180);
close(gcf);
end
