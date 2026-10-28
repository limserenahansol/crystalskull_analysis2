function extract_melittin_lab_meeting(which)
%EXTRACT_MELITTIN_LAB_MEETING  Face-poke-style figures for BL vs injection.
%
% Independent EXTRACT cell IDs (BL and INJ are not the same neurons).
% Injection is cropped to the last N seconds matching baseline duration.
% Traces: EXTRACT /T. Maps: frr /S footprints on contrast-stretched FOV.
%
%   extract_melittin_lab_meeting         % all recordings + cohort PPT figs
%   extract_melittin_lab_meeting(3)     % one job (1-based)
%   extract_melittin_lab_meeting("cohort")

if nargin < 1 || isempty(which)
    which = "all";
end

root = fileparts(mfilename("fullpath"));
addpath(root);

outDir = fullfile(root, "melittin_lab_meeting");
perDir = fullfile(outDir, "per_recording");
avgDir = fullfile(outDir, "cohort_average");
exDir = fullfile(outDir, "example_neurons", "per_recording");
if exist(perDir, "dir") ~= 7
    mkdir(perDir);
end
if exist(avgDir, "dir") ~= 7
    mkdir(avgDir);
end
if exist(exDir, "dir") ~= 7
    mkdir(exDir);
end

jobs = melittin_jobs();
nJobs = size(jobs, 1);

if (ischar(which) || isstring(which)) && strcmpi(strtrim(char(string(which))), "cohort")
    write_cohort_from_saved(outDir, jobs);
    return
end

jobList = 1:nJobs;
if isnumeric(which) && isfinite(which) && which >= 1
    jobList = round(double(which));
    jobList = jobList(jobList >= 1 & jobList <= nJobs);
    if isempty(jobList)
        error("extract_melittin_lab_meeting:BadJob", "Job index out of range 1..%d.", nJobs);
    end
end

nOk = 0;
for ii = 1:numel(jobList)
    j = jobList(ii);
    rec = jobs(j);
    recLabel = rec_label(rec);
    fprintf("\n==== %d/%d  %s ====\n", j, nJobs, recLabel);
    recOut = fullfile(perDir, recLabel);
    exOut = fullfile(exDir, recLabel);
    if exist(recOut, "dir") ~= 7
        mkdir(recOut);
    end
    if exist(exOut, "dir") ~= 7
        mkdir(exOut);
    end
    matFile = fullfile(recOut, "results.mat");
    try
        if exist(matFile, "file") == 2 && exist(fullfile(exOut, "example_traces.png"), "file") == 2 ...
                && exist(fullfile(recOut, "fig4_spatial_response.png"), "file") == 2
            fprintf("  exists; loading %s\n", matFile);
        else
            R = process_one(rec, recOut, exOut);
            save(matFile, "R", "-v7.3");
        end
        nOk = nOk + 1;
    catch ME
        warning("extract_melittin_lab_meeting:JobFailed", "%s: %s", recLabel, ME.message);
        disp(getReport(ME, "extended", "hyperlinks", "off"));
    end
end

if numel(jobList) == nJobs || strcmpi(char(string(which)), "all")
    write_cohort_from_saved(outDir, jobs);
end
fprintf("\nDone  %d recordings processed this call  -> %s\n", nOk, outDir);
end


function jobs = melittin_jobs()
jobs = [ ...
    job_row("melittin_only", "\\\\stg-tnr50a1.stanford.edu\\crystal_skull\\omer_hazon\\meltinonly", "CS1014-1b", "mirror1path", "CS1014");
    job_row("melittin_only", "\\\\stg-tnr50a1.stanford.edu\\crystal_skull\\omer_hazon\\meltinonly", "CS1014-2w", "mirror1path", "CS1014");
    job_row("melittin_only", "\\\\stg-tnr50a1.stanford.edu\\crystal_skull\\omer_hazon\\meltinonly", "CS1014-1b", "mirror2path", "CS1014");
    job_row("melittin_only", "\\\\stg-tnr50a1.stanford.edu\\crystal_skull\\omer_hazon\\meltinonly", "CS1014-2w", "mirror2path", "CS1014");
    job_row("melittin_only", "\\\\stg-tnr50a1.stanford.edu\\crystal_skull\\omer_hazon", "CS0204-1b", "mirror1path", "CS0204");
    job_row("melittin_only", "\\\\stg-tnr50a1.stanford.edu\\crystal_skull\\omer_hazon", "CS0204-2w", "mirror1path", "CS0204");
    job_row("melittin_only", "\\\\stg-tnr50a1.stanford.edu\\crystal_skull\\omer_hazon", "CS0204-1b", "mirror2path", "CS0204");
    job_row("melittin_only", "\\\\stg-tnr50a1.stanford.edu\\crystal_skull\\omer_hazon", "CS0204-2w", "mirror2path", "CS0204");
    job_row("melittin_morphine", "\\\\stg-tnr50a1.stanford.edu\\crystal_skull\\omer_hazon\\morphine_20260408", "1014m01b", "mirror1path", "morphine");
    job_row("melittin_morphine", "\\\\stg-tnr50a1.stanford.edu\\crystal_skull\\omer_hazon\\morphine_20260408", "1014m02w", "mirror1path", "morphine");
    job_row("melittin_morphine", "\\\\stg-tnr50a1.stanford.edu\\crystal_skull\\omer_hazon\\morphine_20260408", "1014m01b", "mirror2path", "morphine");
    job_row("melittin_morphine", "\\\\stg-tnr50a1.stanford.edu\\crystal_skull\\omer_hazon\\morphine_20260408", "1014m02w", "mirror2path", "morphine") ...
    ];
end


function s = job_row(cohort, root, mouse, mirror, subgroup)
s = struct("cohort", char(cohort), "root", char(root), "mouse", char(mouse), ...
    "mirror", char(mirror), "subgroup", char(subgroup));
end


function label = rec_label(rec)
mir = rec.mirror;
if strcmp(mir, "mirror1path")
    mir = "mirror1";
elseif strcmp(mir, "mirror2path")
    mir = "mirror2";
end
label = sprintf("%s_%s_%s", rec.mouse, rec.subgroup, char(mir));
end


function txt = rec_title(rec)
mir = rec.mirror;
if strcmp(mir, "mirror1path")
    mir = "mirror1";
elseif strcmp(mir, "mirror2path")
    mir = "mirror2";
end
txt = sprintf("%s  |  %s  |  %s", rec.mouse, strrep(rec.cohort, "_", " "), char(mir));
end


function R = process_one(rec, recOut, exOut)
blPath = fullfile(rec.root, rec.mouse, "Baseline", rec.mirror);
injPath = fullfile(rec.root, rec.mouse, "Injection", rec.mirror);
Rb = extract_resolve_extract_files(blPath);
Ri = extract_resolve_extract_files(injPath);
if isempty(Rb.frr) || isempty(Ri.frr) || isempty(Rb.summary) || isempty(Ri.summary)
    error("extract_melittin_lab_meeting:MissingFiles", "Missing EXTRACT files for %s.", rec_label(rec));
end

fsBl = double(h5read(Rb.summary, "/framerate"));
fsInj = double(h5read(Ri.summary, "/framerate"));
fsBl = fsBl(1);
fsInj = fsInj(1);

fprintf("  Loading baseline T (chunked 1 s bins)...\n");
Bbl = load_binned_session(Rb.frr, fsBl, inf);
fprintf("  Loading injection T (full)...\n");
Binj = load_binned_session(Ri.frr, fsInj, inf);
nSec = min(Bbl.dur, Binj.dur);
fprintf("  Matching both to the shorter window: %.1f s  (BL %.1f s, INJ %.1f s)\n", ...
    nSec, Bbl.dur, Binj.dur);
Bbl = crop_last_sec(Bbl, nSec);
Binj = crop_last_sec(Binj, nSec);

dffBl = dff_from_bins(Bbl.Tbin);
dffInj = dff_from_bins(Binj.Tbin);
zBl = zscore_cols(dffBl);
zInj = zscore_cols(dffInj);

popDffBl = mean(dffBl, 2, "omitnan");
popDffInj = mean(dffInj, 2, "omitnan");
popFracBl = mean(zBl > 2, 2, "omitnan");
popFracInj = mean(zInj > 2, 2, "omitnan");
meanDffBl = mean(dffBl, 1, "omitnan")';
meanDffInj = mean(dffInj, 1, "omitnan")';

amp = meanDffInj - mean(meanDffBl, "omitnan");
nBins = size(dffInj, 1);
delta = dffInj - mean(meanDffBl, "omitnan");
se = std(delta, 0, 1, "omitnan") ./ sqrt(max(nBins, 1));
tstat = mean(delta, 1, "omitnan") ./ max(se, 1e-8);
pval = 2 * tcdf(-abs(tstat), max(nBins - 1, 1));
pval = pval(:);
isUp = pval < 0.05 & amp > 0;
isDown = pval < 0.05 & amp < 0;

try
    pPop = ranksum(popDffBl(:), popDffInj(:));
catch
    pPop = NaN;
end

flipLR = contains(rec.mirror, "mirror2");
sumimg = [];
try
    if ~isempty(Ri.ds_ext)
        sumimg = h5read(Ri.ds_ext, "/info/summary_image")';
    end
catch
end
if flipLR && ~isempty(sumimg)
    sumimg = fliplr(sumimg);
end

fprintf("  Reconstructing injection footprints...\n");
[Sn, H] = load_frr_spatial(Ri.frr, Binj.nNeurons, 2304);
nPix = size(Sn, 1);
rowS = rem((0:nPix-1)', H) + 1;
colSum = full(sum(Sn, 1));
if flipLR
    cx = full(((H + 1) - rowS)' * Sn) ./ max(colSum, eps);
else
    cx = full((rowS' * Sn) ./ max(colSum, eps));
end
cx = cx(:);
ok = colSum(:) > 0 & isfinite(amp) & isfinite(cx) & isfinite(pval);

ampSig = amp;
ampSig(~(isUp | isDown)) = 0;
respMap = reshape(full(Sn * ampSig), H, H)';
if flipLR
    respMap = fliplr(respMap);
end
showPx = 768;
respSmall = imresize(respMap, [showPx, showPx]);

label = rec_label(rec);
dispLabel = rec_title(rec);
write_rec_figures(recOut, dispLabel, Bbl, Binj, popDffBl, popDffInj, popFracBl, popFracInj, ...
    zInj, amp, isUp, isDown, pPop, respSmall);

examples = pick_examples(amp, isUp, isDown, cx, ok);
for k = 1:numel(examples)
    examples(k).nUp = sum(isUp);
    examples(k).nDown = sum(isDown);
    examples(k).pPop = pPop;
end
write_example_csv(examples, exOut);
traces = load_example_traces(Ri.frr, Binj, [examples.id]);
write_trace_figure(traces, Binj, examples, popDffInj, mean(meanDffBl, "omitnan"), dispLabel, exOut);
write_map_figure(sumimg, Sn, H, amp, isUp, isDown, examples, dispLabel, exOut, flipLR);

R = struct();
R.cohort = rec.cohort;
R.subgroup = rec.subgroup;
R.mouse = rec.mouse;
R.mirror = rec.mirror;
R.label = label;
R.n_neurons_bl = Bbl.nNeurons;
R.n_neurons_inj = Binj.nNeurons;
R.dur_bl = Bbl.dur;
R.dur_inj = Binj.dur;
R.fs_bl = fsBl;
R.fs_inj = fsInj;
R.t_bl = Bbl.t;
R.t_inj = Binj.t;
R.pop_dff_bl = popDffBl;
R.pop_dff_inj = popDffInj;
R.pop_frac_bl = popFracBl;
R.pop_frac_inj = popFracInj;
R.mean_dff_bl = mean(meanDffBl, "omitnan");
R.mean_dff_inj = mean(meanDffInj, "omitnan");
R.d_pop_dff = R.mean_dff_inj - R.mean_dff_bl;
R.frac_up = mean(isUp);
R.frac_down = mean(isDown);
R.n_up = sum(isUp);
R.n_down = sum(isDown);
R.p_pop = pPop;
R.amp = amp(:);
R.pval = pval(:);
R.is_up = isUp(:);
R.is_down = isDown(:);
R.cx = cx(:);
R.mean_dff_inj_cell = meanDffInj(:);
R.resp_map_small = respSmall;
R.sta_up = mean_trace(zInj(:, isUp));
R.sta_down = mean_trace(zInj(:, isDown));
R.sta_all = mean_trace(zInj);
end


function B = load_binned_session(frrFile, fs, useLastSec)
info = h5info(frrFile, "/T");
sz = double(info.Dataspace.Size);
jc = double(h5read(frrFile, "/S/jc"));
nS = numel(jc) - 1;
if sz(1) == nS
    neuDim = 1;
    frDim = 2;
    nNeu = sz(1);
    nFr = sz(2);
elseif sz(2) == nS
    neuDim = 2;
    frDim = 1;
    nNeu = sz(2);
    nFr = sz(1);
else
    neuDim = 1;
    frDim = 2;
    nNeu = sz(1);
    nFr = sz(2);
    warning("extract_melittin_lab_meeting:TShape", ...
        "T size %s does not match S columns %d; assuming neurons x frames.", mat2str(sz), nS);
end

if isfinite(useLastSec) && useLastSec > 0 && useLastSec < (nFr / fs - 0.5)
    take = min(nFr, max(1, round(double(useLastSec) * double(fs))));
    f0 = nFr - take + 1;
else
    f0 = 1;
    take = nFr;
end

binF = max(1, round(double(fs)));
nBins = floor(take / binF);
if nBins < 5
    error("extract_melittin_lab_meeting:TooShort", "Only %d 1-s bins after crop.", nBins);
end

Tbin = zeros(nBins, nNeu, "single");
chunkBins = 40;
b = 1;
while b <= nBins
    nTake = min(chunkBins, nBins - b + 1);
    nFrTake = nTake * binF;
    fr = f0 + (b - 1) * binF;
    if neuDim == 1
        chunk = double(h5read(frrFile, "/T", [1, fr], [nNeu, nFrTake]));
        for k = 1:nTake
            sl = (k - 1) * binF + 1 : k * binF;
            Tbin(b + k - 1, :) = mean(chunk(:, sl), 2, "omitnan")';
        end
    else
        chunk = double(h5read(frrFile, "/T", [fr, 1], [nFrTake, nNeu]));
        for k = 1:nTake
            sl = (k - 1) * binF + 1 : k * binF;
            Tbin(b + k - 1, :) = mean(chunk(sl, :), 1, "omitnan");
        end
    end
    b = b + nTake;
end

B.Tbin = Tbin;
B.nNeurons = nNeu;
B.nFramesFile = nFr;
B.nFramesUsed = nBins * binF;
B.f0 = f0;
B.take = take;
B.binF = binF;
B.fs = double(fs);
B.dur = B.nFramesUsed / B.fs;
B.t = ((0:nBins-1)' + 0.5);
B.neuDim = neuDim;
B.frrFile = frrFile;
end


function B = crop_last_sec(B, nSec)
nKeep = min(size(B.Tbin, 1), max(5, round(double(nSec))));
nDrop = size(B.Tbin, 1) - nKeep;
if nDrop > 0
    B.Tbin = B.Tbin(end-nKeep+1:end, :);
    framesDrop = nDrop * B.binF;
    B.f0 = B.f0 + framesDrop;
end
B.take = nKeep * B.binF;
B.nFramesUsed = B.take;
B.dur = nKeep;
B.t = ((0:nKeep-1)' + 0.5);
end


function dff = dff_from_bins(Tbin)
F0 = prctile(double(Tbin), 8, 1);
F0 = max(F0, 0.5);
dff = (double(Tbin) - F0) ./ F0;
end


function z = zscore_cols(x)
z = (x - mean(x, 1, "omitnan")) ./ max(std(x, 0, 1, "omitnan"), 1e-8);
end


function mu = mean_trace(z)
if isempty(z)
    mu = [];
    return
end
mu = mean(z, 2, "omitnan");
end


function [Sn, H] = load_frr_spatial(frrFile, nNeurons, framePx)
ir = double(h5read(frrFile, "/S/ir"));
jc = double(h5read(frrFile, "/S/jc"));
data = double(h5read(frrFile, "/S/data"));
nCol = numel(jc) - 1;
nPix = framePx * framePx;
colIdx = zeros(jc(end), 1);
for c = 1:nCol
    colIdx(jc(c)+1:jc(c+1)) = c;
end
Ss = sparse(ir + 1, colIdx, data, nPix, nCol);
if nCol > nNeurons
    Ss = Ss(:, 1:nNeurons);
elseif nCol < nNeurons
    Ss(nPix, nNeurons) = 0;
end
cs = full(sum(Ss, 1));
cs(cs == 0) = 1;
Sn = Ss ./ cs;
H = framePx;
end


function traces = load_example_traces(frrFile, B, ids)
nEx = numel(ids);
traces = zeros(B.take, nEx);
for k = 1:nEx
    id = ids(k);
    if B.neuDim == 1
        tr = double(h5read(frrFile, "/T", [id, B.f0], [1, B.take]));
        traces(:, k) = tr(:);
    else
        tr = double(h5read(frrFile, "/T", [B.f0, id], [B.take, 1]));
        traces(:, k) = tr(:);
    end
end
end


function examples = pick_examples(amp, isUp, isDown, cx, ok)
examples = struct("id", {}, "kind", {}, "tag", {}, "amp", {}, "cx", {}, "color", {}, ...
    "nUp", {}, "nDown", {}, "pPop", {});
colUp = [0.85, 0.20, 0.15];
colDn = [0.10, 0.35, 0.85];
upMask = isUp & ok;
dnMask = isDown & ok;
if sum(upMask) < 2
    upMask = amp > 0 & ok;
end
if sum(dnMask) < 2
    dnMask = amp < 0 & ok;
end
examples = add_pair(examples, amp, upMask, cx, "up", colUp, true);
examples = add_pair(examples, amp, dnMask, cx, "down", colDn, false);
end


function examples = add_pair(examples, amp, mask, cx, kind, color, wantPositive)
ids = find(mask);
if numel(ids) < 2
    ids = find(isfinite(amp) & isfinite(cx));
end
if numel(ids) < 2
    error("extract_melittin_lab_meeting:TooFew", "Need at least 2 %s cells.", kind);
end
if wantPositive
    score = amp;
else
    score = -amp;
end
medX = median(cx(ids));
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
if strcmp(kind, "up")
    tagAnt = "U1 anterior";
    tagPost = "U2 posterior";
else
    tagAnt = "D1 anterior";
    tagPost = "D2 posterior";
end
examples(end+1) = pack_example(idAnt, kind, tagAnt, amp(idAnt), cx(idAnt), color);
examples(end+1) = pack_example(idPost, kind, tagPost, amp(idPost), cx(idPost), color);
end


function s = pack_example(id, kind, tag, amp, cx, color)
s = struct("id", id, "kind", kind, "tag", tag, "amp", amp, "cx", cx, "color", color, ...
    "nUp", 0, "nDown", 0, "pPop", NaN);
end


function write_example_csv(examples, outDir)
fid = fopen(fullfile(outDir, "example_neurons.csv"), "w");
fprintf(fid, "tag,neuron_id,kind,amp,centroid_x\n");
for i = 1:numel(examples)
    e = examples(i);
    fprintf(fid, "%s,%d,%s,%.5f,%.2f\n", char(e.tag), e.id, e.kind, e.amp, e.cx);
end
fclose(fid);
end


function write_rec_figures(recOut, dispLabel, Bbl, Binj, popDffBl, popDffInj, popFracBl, popFracInj, ...
    zInj, amp, isUp, isDown, pPop, respSmall)
colBl = [0.20, 0.45, 0.75];
colInj = [0.90, 0.40, 0.12];
colUp = [0.85, 0.20, 0.15];
colDn = [0.10, 0.35, 0.85];

figure("Position", [20, 20, 1500, 560], "Color", "w", "Visible", "off");
sgtitle(sprintf("Fig1  |  Population activity  |  %s", dispLabel), ...
    "FontSize", 13, "FontWeight", "bold", "Interpreter", "none");
subplot(1, 2, 1);
plot(Bbl.t, popDffBl, "Color", colBl, "LineWidth", 1.4); hold on;
plot(Binj.t, popDffInj, "Color", colInj, "LineWidth", 1.4);
xlabel("Time in matched window (s)");
ylabel("Population mean dF/F");
title("Blue = baseline  |  orange = injection (matched to shorter duration)");
box off;
legend({"Baseline", "Injection"}, "Location", "best");
subplot(1, 2, 2);
plot(Bbl.t, popFracBl, "Color", colBl, "LineWidth", 1.4); hold on;
plot(Binj.t, popFracInj, "Color", colInj, "LineWidth", 1.4);
xlabel("Time in matched window (s)");
ylabel("Fraction active (z>2)");
title(sprintf("Population p = %.3g (ranksum of mean dF/F bins)", pPop));
box off;
exportgraphics(gcf, fullfile(recOut, "fig1_population_psth.png"), "Resolution", 180);
close(gcf);

figure("Position", [20, 20, 900, 520], "Color", "w", "Visible", "off");
tInj = Binj.t;
yu = mean_trace(zInj(:, isUp));
yd = mean_trace(zInj(:, isDown));
ya = mean_trace(zInj);
if isempty(yu)
    yu = nan(size(tInj));
end
if isempty(yd)
    yd = nan(size(tInj));
end
plot(tInj, yu, "Color", colUp, "LineWidth", 2.2); hold on;
plot(tInj, yd, "Color", colDn, "LineWidth", 2.2);
plot(tInj, ya, "k", "LineWidth", 1.6);
yline(0, "k:", "LineWidth", 0.8);
xlabel("Time in injection window (s)");
ylabel("Mean z (injection cells)");
title(sprintf("Fig2  |  Up / down / all injection cells  |  %s", dispLabel), "Interpreter", "none");
legend({"Up vs BL mean", "Down vs BL mean", "All injection cells"}, "Location", "best");
box off;
exportgraphics(gcf, fullfile(recOut, "fig2_neuron_sta.png"), "Resolution", 180);
close(gcf);

figure("Position", [20, 20, 1400, 500], "Color", "w", "Visible", "off");
sgtitle(sprintf("Fig3  |  Injection cells vs baseline session mean  |  %s", dispLabel), ...
    "FontSize", 13, "FontWeight", "bold", "Interpreter", "none");
subplot(1, 3, 1);
b = bar([1 2], [mean(isUp), mean(isDown)] * 100, 0.55, "EdgeColor", "k");
b.FaceColor = "flat";
b.CData = [colUp; colDn];
set(gca, "XTick", 1:2, "XTickLabel", {"Up", "Down"});
ylabel("% of injection neurons");
title(sprintf("p<0.05 vs BL mean dF/F"));
box off;
subplot(1, 3, 2);
histogram(amp, 40, "FaceColor", [0.45, 0.45, 0.55], "EdgeColor", "none");
xline(0, "k-", "LineWidth", 1.2);
xlabel("amp = mean dF/F_inj - mean(BL cells)");
ylabel("# injection neurons");
title("All injection cells");
box off;
subplot(1, 3, 3);
text(0.08, 0.72, sprintf("Up: %d (%.1f%%)\nDown: %d (%.1f%%)\nNeurons inj: %d\nNeurons BL: %d\nPop p = %.3g", ...
    sum(isUp), 100*mean(isUp), sum(isDown), 100*mean(isDown), size(zInj, 2), Bbl.nNeurons, pPop), ...
    "FontSize", 13, "Interpreter", "none");
axis off;
exportgraphics(gcf, fullfile(recOut, "fig3_responsive_cells.png"), "Resolution", 180);
close(gcf);

halfN = 128;
cmapRb = [[linspace(0, 1, halfN), ones(1, halfN)]' ...
    [linspace(1, 1, halfN), linspace(1, 0, halfN)]' ...
    [ones(1, halfN), linspace(1, 0, halfN)]'];
figure("Position", [20, 20, 820, 760], "Color", "w", "Visible", "off");
sigAbs = abs(respSmall(respSmall ~= 0));
climv = prctile(sigAbs, 95);
if isempty(climv) || ~isfinite(climv) || climv == 0
    climv = 1;
end
imagesc(respSmall);
colormap(gca, cmapRb);
clim([-climv, climv]);
axis image off;
add_ap_labels(gca);
title(sprintf("Fig4  |  Injection FOV  red=up  blue=down  |  %s", dispLabel), ...
    "Interpreter", "none", "FontSize", 12);
cb = colorbar;
cb.Label.String = "amp (inj mean dF/F - BL session)";
exportgraphics(gcf, fullfile(recOut, "fig4_spatial_response.png"), "Resolution", 180);
close(gcf);

figure("Position", [20, 20, 900, 420], "Color", "w", "Visible", "off");
sgtitle(sprintf("Fig5  |  Session length and cell counts  |  %s", dispLabel), ...
    "FontSize", 13, "FontWeight", "bold", "Interpreter", "none");
subplot(1, 2, 1);
bar([1 2], [Bbl.dur, Binj.dur], 0.55, "FaceColor", [0.4, 0.45, 0.75], "EdgeColor", "none");
set(gca, "XTick", 1:2, "XTickLabel", {"Baseline", "Injection"});
ylabel("Duration used (s)");
title("Matched injection tail");
box off;
subplot(1, 2, 2);
bar([1 2], [Bbl.nNeurons, Binj.nNeurons], 0.55, "FaceColor", [0.4, 0.4, 0.55], "EdgeColor", "none");
set(gca, "XTick", 1:2, "XTickLabel", {"Baseline", "Injection"});
ylabel("# EXTRACT neurons");
title("Independent cell IDs (not matched)");
box off;
exportgraphics(gcf, fullfile(recOut, "fig5_session_stats.png"), "Resolution", 180);
close(gcf);
end


function write_trace_figure(traces, Binj, examples, popDffInj, blMean, label, outDir)
nEx = numel(examples);
t = (0:size(traces, 1)-1)' / Binj.fs;
winSec = 45;
binF = Binj.binF;
nB = floor(size(traces, 1) / binF);
proxy = zeros(nB, 1);
for b = 1:nB
    sl = (b - 1) * binF + 1 : b * binF;
    chunk = traces(sl, :);
    proxy(b) = mean(abs(chunk - median(chunk, 1, "omitnan")), "all", "omitnan");
end
tStart = choose_high_activity_window(Binj.t(1:min(nB, numel(Binj.t))), proxy, winSec, t(end));
tEnd = tStart + winSec;
inWin = t >= tStart & t < tEnd;

figure("Position", [20, 20, 1760, 980], "Color", "w", "Visible", "off");
sgtitle({sprintf("Example calcium traces   |   %s", label), ...
    "Injection EXTRACT traces (left) and 1 s dF/F of the same cell (right). Cells are injection-only (not paired to baseline)."}, ...
    "FontSize", 13, "FontWeight", "bold", "Interpreter", "none");

for i = 1:nEx
    e = examples(i);
    tr = traces(inWin, i);
    tr = tr - median(tr, "omitnan");
    ax = subplot(nEx, 2, 2*i-1);
    plot(ax, t(inWin), tr, "Color", e.color, "LineWidth", 1.15);
    ylabel(ax, "dF/F-like");
    title(ax, sprintf("%s   |   neuron %d   |   amp = %+.3f", char(e.tag), e.id, e.amp), ...
        "Interpreter", "none", "Color", e.color, "FontSize", 11);
    box(ax, "off");
    xlim(ax, [tStart, tEnd]);
    if i == nEx
        xlabel(ax, "Time in injection window (s)");
    else
        set(ax, "XTickLabel", []);
    end

    ax2 = subplot(nEx, 2, 2*i);
    hold(ax2, "on");
    yline(ax2, 0, "k:", "LineWidth", 0.8);
    yline(ax2, blMean, "--", "Color", [0.20, 0.45, 0.75], "LineWidth", 1.2);
    % 1 s bins of this cell: reconstruct from high-res by binning
    binF = Binj.binF;
    nB = floor(size(traces, 1) / binF);
    tb = ((0:nB-1)' + 0.5);
    F0 = max(prctile(traces(:, i), 8), 0.5);
    dff1 = zeros(nB, 1);
    for b = 1:nB
        sl = (b-1)*binF+1 : b*binF;
        dff1(b) = mean(traces(sl, i), "omitnan") / F0 - 1;
    end
    plot(ax2, tb, dff1, "Color", e.color, "LineWidth", 1.6);
    ylabel(ax2, "dF/F");
    title(ax2, "Injection 1 s bins   |   dashed = BL session mean", "FontSize", 10);
    box(ax2, "off");
    xlim(ax2, [tb(1), tb(end)]);
    if i == nEx
        xlabel(ax2, "Time in injection window (s)");
    else
        set(ax2, "XTickLabel", []);
    end
end
exportgraphics(gcf, fullfile(outDir, "example_traces.png"), "Resolution", 200);
close(gcf);
end


function tStart = choose_high_activity_window(tBin, popDff, winSec, tMax)
best = 0;
bestMu = -inf;
for i = 1:numel(tBin)
    t0 = max(tBin(i) - 2, 0);
    if t0 + winSec > tMax
        break
    end
    mu = mean(popDff(tBin >= t0 & tBin < t0 + winSec), "omitnan");
    if mu > bestMu
        bestMu = mu;
        best = t0;
    end
end
tStart = best;
end


function write_map_figure(sumimg, Sn, H, amp, isUp, isDown, examples, label, outDir, flipLR)
showPx = 768;
halfN = 128;
cmapRb = [[linspace(0, 1, halfN), ones(1, halfN)]' ...
    [linspace(0, 1, halfN), linspace(1, 0, halfN)]' ...
    [ones(1, halfN), linspace(1, 0, halfN)]'];

ampSig = amp;
ampSig(~(isUp | isDown)) = 0;
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
    "Injection crystal-skull FOV. P = posterior, A = anterior. Red = up vs BL mean, blue = down."}, ...
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
title(ax2, sprintf("Significant injection cells   up=%d   down=%d", sum(isUp), sum(isDown)), ...
    "FontSize", 13);
add_ap_labels(ax2);
cb = colorbar(ax2);
cb.Label.String = "amp (inj - BL session, dF/F)";

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
try
    img = adapthisteq(img, "NumTiles", [16, 16], "ClipLimit", 0.04, "Range", "original");
catch
end
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
        try
            [B, ~] = bwboundaries(mask, "noholes");
            for b = 1:numel(B)
                plot(ax, B{b}(:, 2), B{b}(:, 1), "Color", e.color, "LineWidth", 1.8);
            end
        catch
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


function write_cohort_from_saved(outDir, jobs)
perDir = fullfile(outDir, "per_recording");
avgDir = fullfile(outDir, "cohort_average");
exDir = fullfile(outDir, "example_neurons", "per_recording");
results = {};
labels = {};
cohorts = {};
for j = 1:size(jobs, 1)
    rec = jobs(j);
    recLabel = rec_label(rec);
    matFile = fullfile(perDir, recLabel, "results.mat");
    if exist(matFile, "file") ~= 2
        continue
    end
    S = load(matFile, "R");
    results{end+1} = S.R; %#ok<AGROW>
    labels{end+1} = recLabel; %#ok<AGROW>
    cohorts{end+1} = rec.cohort; %#ok<AGROW>
end
if isempty(results)
    warning("extract_melittin_lab_meeting:NoResults", "No saved results yet; skip cohort figures.");
    return
end
write_index_csv(outDir, results);
write_cohort_five(results, avgDir, "all", "All recordings");
only = filter_cohort(results, "melittin_only");
morph = filter_cohort(results, "melittin_morphine");
if ~isempty(only)
    write_cohort_five(only, avgDir, "melittin_only", "Melittin only  (CS1014 + CS0204)");
end
if ~isempty(morph)
    write_cohort_five(morph, avgDir, "morphine", "Melittin + morphine");
end
write_takehome(results, avgDir);
write_overview(exDir, labels, cohorts, fullfile(outDir, "example_neurons"));
write_master_example_csv(exDir, labels, cohorts, outDir);
write_readme(outDir);
end


function out = filter_cohort(results, name)
out = {};
for i = 1:numel(results)
    if strcmp(char(string(results{i}.cohort)), name)
        out{end+1} = results{i}; %#ok<AGROW>
    end
end
end


function write_cohort_five(results, avg, tag, titleStem)
if nargin < 3 || isempty(tag)
    tag = "all";
end
if nargin < 4 || isempty(titleStem)
    titleStem = "All recordings";
end
tag = char(string(tag));
if strcmp(tag, "all")
    suffix = "";
else
    suffix = sprintf("_%s", tag);
end
n = numel(results);
nq = 400;
tq = linspace(0, 1, nq);
psthBl = nan(n, nq);
psthInj = nan(n, nq);
fracBl = nan(n, nq);
fracInj = nan(n, nq);
staU = nan(n, nq);
staD = nan(n, nq);
staA = nan(n, nq);
fe = nan(n, 1);
fi = nan(n, 1);
pp = nan(n, 1);
nn = nan(n, 1);
dPop = nan(n, 1);
colUp = [0.85, 0.20, 0.15];
colDn = [0.10, 0.35, 0.85];
colBl = [0.20, 0.45, 0.75];
colInj = [0.90, 0.40, 0.12];

for i = 1:n
    R = results{i};
    tb = (R.t_bl - R.t_bl(1)) / max(R.t_bl(end) - R.t_bl(1), eps);
    ti = (R.t_inj - R.t_inj(1)) / max(R.t_inj(end) - R.t_inj(1), eps);
    psthBl(i, :) = interp1(tb, R.pop_dff_bl, tq, "linear", NaN);
    psthInj(i, :) = interp1(ti, R.pop_dff_inj, tq, "linear", NaN);
    fracBl(i, :) = interp1(tb, R.pop_frac_bl, tq, "linear", NaN);
    fracInj(i, :) = interp1(ti, R.pop_frac_inj, tq, "linear", NaN);
    if ~isempty(R.sta_up)
        tu = linspace(0, 1, numel(R.sta_up));
        staU(i, :) = interp1(tu, R.sta_up, tq, "linear", NaN);
    end
    if ~isempty(R.sta_down)
        tu = linspace(0, 1, numel(R.sta_down));
        staD(i, :) = interp1(tu, R.sta_down, tq, "linear", NaN);
    end
    if ~isempty(R.sta_all)
        tu = linspace(0, 1, numel(R.sta_all));
        staA(i, :) = interp1(tu, R.sta_all, tq, "linear", NaN);
    end
    fe(i) = R.frac_up;
    fi(i) = R.frac_down;
    pp(i) = R.p_pop;
    nn(i) = R.n_neurons_inj;
    dPop(i) = R.d_pop_dff;
end
sem = @(x) std(x, 0, 1, "omitnan") ./ sqrt(max(sum(isfinite(x), 1), 1));

figure("Position", [30, 30, 1500, 560], "Color", "w", "Visible", "off");
sgtitle(sprintf("Fig1  |  %s (n=%d): population BL vs injection", titleStem, n), ...
    "FontSize", 13, "FontWeight", "bold");
subplot(1, 2, 1);
plot(tq, psthBl', "Color", [0.75, 0.80, 0.90], "LineWidth", 0.7); hold on;
plot(tq, psthInj', "Color", [0.93, 0.80, 0.72], "LineWidth", 0.7);
muB = mean(psthBl, 1, "omitnan"); seB = sem(psthBl);
muI = mean(psthInj, 1, "omitnan"); seI = sem(psthInj);
fill([tq fliplr(tq)], [muB-seB fliplr(muB+seB)], colBl, "FaceAlpha", 0.25, "EdgeColor", "none");
fill([tq fliplr(tq)], [muI-seI fliplr(muI+seI)], colInj, "FaceAlpha", 0.25, "EdgeColor", "none");
plot(tq, muB, "Color", colBl, "LineWidth", 2.5);
plot(tq, muI, "Color", colInj, "LineWidth", 2.5);
xlabel("Normalized time in matched window");
ylabel("Population mean dF/F");
title("Gray-tint = each recording  |  bold = mean \pm SEM");
box off;
hBl = plot(nan, nan, "Color", colBl, "LineWidth", 2.5);
hInj = plot(nan, nan, "Color", colInj, "LineWidth", 2.5);
legend([hBl, hInj], {"Baseline mean", "Injection mean"}, "Location", "best");
subplot(1, 2, 2);
plot(tq, fracBl', "Color", [0.75, 0.80, 0.90], "LineWidth", 0.7); hold on;
plot(tq, fracInj', "Color", [0.93, 0.80, 0.72], "LineWidth", 0.7);
muB = mean(fracBl, 1, "omitnan"); seB = sem(fracBl);
muI = mean(fracInj, 1, "omitnan"); seI = sem(fracInj);
fill([tq fliplr(tq)], [muB-seB fliplr(muB+seB)], colBl, "FaceAlpha", 0.25, "EdgeColor", "none");
fill([tq fliplr(tq)], [muI-seI fliplr(muI+seI)], colInj, "FaceAlpha", 0.25, "EdgeColor", "none");
plot(tq, muB, "Color", colBl, "LineWidth", 2.5);
plot(tq, muI, "Color", colInj, "LineWidth", 2.5);
xlabel("Normalized time in matched window");
ylabel("Fraction active (z>2)");
title("Population recruitment");
box off;
exportgraphics(gcf, fullfile(avg, sprintf("fig1_population_psth%s.png", suffix)), "Resolution", 200);
close(gcf);

figure("Position", [30, 30, 900, 520], "Color", "w", "Visible", "off");
plot(tq, mean(staU, 1, "omitnan"), "Color", colUp, "LineWidth", 2.2); hold on;
plot(tq, mean(staD, 1, "omitnan"), "Color", colDn, "LineWidth", 2.2);
plot(tq, mean(staA, 1, "omitnan"), "k", "LineWidth", 1.8);
yline(0, "k:", "LineWidth", 0.8);
xlabel("Normalized time in injection window");
ylabel("Mean z (injection cells)");
title(sprintf("Fig2  |  Cohort-mean up / down / all  |  %s (n=%d)", titleStem, n));
legend({"Up cells (mean of movies)", "Down cells", "All injection cells"}, "Location", "best");
box off;
exportgraphics(gcf, fullfile(avg, sprintf("fig2_neuron_sta%s.png", suffix)), "Resolution", 200);
close(gcf);

figure("Position", [30, 30, 1400, 500], "Color", "w", "Visible", "off");
sgtitle(sprintf("Fig3  |  %s  |  injection cells vs baseline session mean", titleStem), "FontSize", 13, "FontWeight", "bold");
subplot(1, 3, 1);
b = bar([1 2], [mean(fe), mean(fi)] * 100, 0.55, "EdgeColor", "k"); hold on;
b.FaceColor = "flat";
b.CData = [colUp; colDn];
errorbar([1 2], [mean(fe), mean(fi)] * 100, [std(fe), std(fi)] * 100, "k.", "LineWidth", 1.4);
set(gca, "XTick", 1:2, "XTickLabel", {"Up", "Down"});
ylabel("% of injection neurons (mean \pm SD across recordings)");
title(sprintf("n = %d recordings", n));
box off;
subplot(1, 3, 2);
plot(100*fe, 100*fi, "ko", "MarkerFaceColor", [0.3, 0.3, 0.3]); hold on;
xlabel("% up");
ylabel("% down");
axis equal;
title("Each recording");
box off;
grid on;
subplot(1, 3, 3);
histogram(log10(max(pp, realmin)), "FaceColor", [0.5, 0.5, 0.7], "EdgeColor", "none");
xlabel("log10 population p (BL vs INJ mean dF/F)");
ylabel("# recordings");
title(sprintf("%d / %d recordings p<0.05", sum(pp < 0.05), n));
box off;
exportgraphics(gcf, fullfile(avg, sprintf("fig3_responsive_cells%s.png", suffix)), "Resolution", 200);
close(gcf);

figure("Position", [20, 20, 1600, 900], "Color", "w", "Visible", "off");
sgtitle({sprintf("Fig4  |  %s  |  injection FOV  (red = up vs BL mean, blue = down)", titleStem), ...
    "Each panel is one recording; FOVs are not spatially aligned. Independent EXTRACT cells."}, ...
    "FontSize", 12, "FontWeight", "bold");
ncol = 4;
nrow = ceil(n / ncol);
halfN = 128;
cmapRb = [[linspace(0, 1, halfN), ones(1, halfN)]' ...
    [linspace(0, 1, halfN), linspace(1, 0, halfN)]' ...
    [ones(1, halfN), linspace(1, 0, halfN)]'];
for i = 1:n
    subplot(nrow, ncol, i);
    R = results{i};
    M = R.resp_map_small;
    climv = prctile(abs(M(M ~= 0)), 98);
    if isempty(climv) || ~isfinite(climv) || climv == 0
        climv = 1;
    end
    imagesc(M);
    colormap(gca, cmapRb);
    clim([-climv, climv]);
    axis image off;
    title(strrep(R.label, "_", "\_"), "FontSize", 8, "Interpreter", "tex");
end
exportgraphics(gcf, fullfile(avg, sprintf("fig4_spatial_response%s.png", suffix)), "Resolution", 180);
close(gcf);

figure("Position", [30, 30, 1400, 500], "Color", "w", "Visible", "off");
sgtitle(sprintf("Fig5  |  %s  |  cell counts and population \\Delta dF/F", titleStem), "FontSize", 13, "FontWeight", "bold");
subplot(1, 2, 1);
bar(1:n, nn, 0.7, "FaceColor", [0.4, 0.45, 0.75], "EdgeColor", "none");
ylabel("# injection EXTRACT neurons");
xlabel("Recording");
title(sprintf("Neurons per injection movie (mean %.0f)", mean(nn)));
box off;
subplot(1, 2, 2);
bar(1:n, dPop, 0.7, "FaceColor", [0.85, 0.45, 0.20], "EdgeColor", "none");
yline(0, "k-", "LineWidth", 1);
ylabel("\Delta population mean dF/F (inj - BL)");
xlabel("Recording");
title("Session-level change");
box off;
exportgraphics(gcf, fullfile(avg, sprintf("fig5_session_stats%s.png", suffix)), "Resolution", 200);
close(gcf);
end


function write_takehome(results, avg)
n = numel(results);
d = nan(n, 1);
grp = strings(n, 1);
mouse = strings(n, 1);
mir = strings(n, 1);
for i = 1:n
    R = results{i};
    d(i) = R.d_pop_dff;
    grp(i) = string(R.cohort);
    mouse(i) = string(R.mouse);
    mir(i) = string(R.mirror);
end

figure("Position", [30, 30, 1200, 560], "Color", "w", "Visible", "off");
sgtitle("Take-home  |  Injection - baseline population mean dF/F", "FontSize", 14, "FontWeight", "bold");
hold on;
x = 1:n;
colA = [0.13, 0.55, 0.70];
colB = [0.90, 0.40, 0.12];
for i = 1:n
    if grp(i) == "melittin_morphine"
        c = colB;
    else
        c = colA;
    end
    bar(x(i), d(i), 0.7, "FaceColor", c, "EdgeColor", "none");
end
yline(0, "k-", "LineWidth", 1);
labs = cellfun(@(R) strrep(char(string(R.label)), "_", " "), results, "UniformOutput", false);
set(gca, "XTick", 1:n, "XTickLabel", labs, "XTickLabelRotation", 35);
ylabel("\Delta pop. mean dF/F");
hA = plot(nan, nan, "s", "MarkerFaceColor", colA, "MarkerEdgeColor", colA);
hB = plot(nan, nan, "s", "MarkerFaceColor", colB, "MarkerEdgeColor", colB);
legend([hA, hB], {"melittin only", "melittin + morphine"}, "Location", "best");
box off;
exportgraphics(gcf, fullfile(avg, "fig_takehome_delta.png"), "Resolution", 200);
close(gcf);

% Paired CS1014 vs morphine (same biological mouse + mirror when names match 1b/2w)
pairs = {
    "CS1014-1b", "1014m01b";
    "CS1014-2w", "1014m02w"
    };
figure("Position", [30, 30, 900, 520], "Color", "w", "Visible", "off");
hold on;
nLine = 0;
for p = 1:size(pairs, 1)
    for k = 1:n
        if ~strcmp(char(mouse(k)), char(string(pairs{p, 1})))
            continue
        end
        hit = 0;
        for j = 1:n
            if strcmp(char(mouse(j)), char(string(pairs{p, 2}))) && strcmp(char(mir(j)), char(mir(k)))
                hit = j;
                break
            end
        end
        if hit == 0
            continue
        end
        nLine = nLine + 1;
        plot([1 2], [d(k), d(hit)], "-o", "Color", [0.45, 0.45, 0.45], "LineWidth", 1.3, ...
            "MarkerFaceColor", [0.2, 0.2, 0.2]);
        text(2.05, d(hit), sprintf("%s %s", char(pairs{p, 1}), short_mirror(mir(k))), ...
            "FontSize", 9, "Interpreter", "none");
    end
end
xlim([0.5, 2.8]);
set(gca, "XTick", [1 2], "XTickLabel", {"melittin only", "melittin + morphine"});
ylabel("\Delta pop. mean dF/F (inj - BL)");
title(sprintf("Paired CS1014 mice (n=%d mouse x mirror)", nLine));
yline(0, "k:", "LineWidth", 0.8);
box off;
exportgraphics(gcf, fullfile(avg, "fig_takehome_paired.png"), "Resolution", 200);
close(gcf);
end


function s = short_mirror(mir)
s = char(mir);
s = strrep(s, "path", "");
end


function write_overview(exDir, allLabels, allCohorts, outDir)
n = numel(allLabels);
if n == 0
    return
end
ncol = 4;
nrow = ceil(n / ncol);
figure("Position", [20, 20, 2000, 1100], "Color", "w", "Visible", "off");
sgtitle("Example neurons on each crystal-skull FOV   |   injection movies, BL vs INJ unmatched cells", ...
    "FontSize", 16, "FontWeight", "bold");
for i = 1:n
    recLabel = allLabels{i};
    fovFile = fullfile(exDir, recLabel, "example_fov.png");
    ax = subplot(nrow, ncol, i);
    if exist(fovFile, "file") == 2
        imshow(imread(fovFile));
    else
        axis(ax, "off");
    end
    title(ax, sprintf("%s  %s", allCohorts{i}, strrep(recLabel, "_", " ")), ...
        "FontSize", 8, "Interpreter", "none");
end
exportgraphics(gcf, fullfile(outDir, "example_fov_overview.png"), "Resolution", 180);
close(gcf);
end


function write_master_example_csv(exDir, labels, cohorts, outDir)
fid = fopen(fullfile(outDir, "example_neurons", "all_example_neurons.csv"), "w");
fprintf(fid, "cohort,label,tag,neuron_id,kind,amp,centroid_x\n");
for i = 1:numel(labels)
    csvFile = fullfile(exDir, labels{i}, "example_neurons.csv");
    if exist(csvFile, "file") ~= 2
        continue
    end
    T = readtable(csvFile, "Delimiter", ",");
    for k = 1:height(T)
        fprintf(fid, "%s,%s,%s,%d,%s,%.5f,%.2f\n", cohorts{i}, labels{i}, ...
            char(string(T.tag{k})), T.neuron_id(k), char(string(T.kind{k})), T.amp(k), T.centroid_x(k));
    end
end
fclose(fid);
end


function write_index_csv(outDir, results)
fid = fopen(fullfile(outDir, "recording_index.csv"), "w");
fprintf(fid, ['cohort,mouse,mirror,label,n_neurons_bl,n_neurons_inj,dur_bl,dur_inj,' ...
    'fs_bl,fs_inj,frac_up,frac_down,p_pop,d_pop_dff,mean_dff_bl,mean_dff_inj,n_up,n_down\n']);
for i = 1:numel(results)
    R = results{i};
    fprintf(fid, "%s,%s,%s,%s,%d,%d,%.3f,%.3f,%.4f,%.4f,%.5f,%.5f,%.3e,%.5f,%.5f,%.5f,%d,%d\n", ...
        R.cohort, R.mouse, R.mirror, R.label, R.n_neurons_bl, R.n_neurons_inj, ...
        R.dur_bl, R.dur_inj, R.fs_bl, R.fs_inj, R.frac_up, R.frac_down, R.p_pop, ...
        R.d_pop_dff, R.mean_dff_bl, R.mean_dff_inj, R.n_up, R.n_down);
end
fclose(fid);
end


function write_readme(outDir)
fid = fopen(fullfile(outDir, "README.txt"), "w");
fprintf(fid, ['Melittin lab-meeting figure pack (face-poke visual analog)\n' ...
    '==============================================================\n\n' ...
    'This is NOT poke-aligned. Baseline and injection are independently extracted.\n' ...
    'Do not treat injection neuron #k as the same cell as baseline neuron #k.\n\n' ...
    'Injection window: last N seconds of BOTH arms, N = min(baseline, injection) duration.\n' ...
    'Per injection neuron: amp = mean dF/F - mean of all baseline neurons mean dF/F.\n' ...
    'Up / down: two-sided t-test of 1 s dF/F bins vs that BL session mean, p<0.05.\n' ...
    'Population p: ranksum of 1 s population mean dF/F, BL vs INJ.\n' ...
    'mirror2 is flipped so P is left and A is right, same as face-poke maps.\n']);
fclose(fid);
end
