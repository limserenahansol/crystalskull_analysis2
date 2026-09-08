function extract_stim_write_spatial_maps(outRoot, extRoot, stimMatFile, params)
%EXTRACT_STIM_WRITE_SPATIAL_MAPS  Rebuild fig4 as a real poke response map.
%
%   extract_stim_write_spatial_maps(outRoot, extRoot, stimMatFile, params)
%
% For each movie folder under outRoot, reload EXTRACT traces + spatial
% footprints, recompute peri-poke amplitude with the same windows as the
% original run, and write fig4_spatial_response.png (anatomy | all cells |
% significant cells). Also copies maps into lab_meeting/per_recording and
% writes a cohort gallery.
%
% params.event_field  'peak' | 'pre' | 'post'
% params.peri_win, base_win, resp_win, alpha, fps_movie, frame_px

if nargin < 4 || isempty(params)
    params = struct();
end
eventField = getp(params, "event_field", "peak");
periWin = getp(params, "peri_win", [-2, 3]);
baseWin = getp(params, "base_win", [-2, 0]);
respWin = getp(params, "resp_win", [0.4, 1.2]);
alpha = getp(params, "alpha", 0.05);
fpsMovie = getp(params, "fps_movie", 33);
framePx = getp(params, "frame_px", 2304);
onlyMovies = getp(params, "only_movies", {});

St = load(stimMatFile);
fnames = cellfun(@char, St.fnames, "UniformOutput", false);
switch lower(eventField)
    case "peak"
        anchors = St.peak_idx;
        efLbl = "peak_idx";
    case "pre"
        anchors = St.peak_pre_idx;
        efLbl = "peak_pre_idx";
    case "post"
        anchors = St.peak_post_idx;
        efLbl = "peak_post_idx";
    otherwise
        error("extract_stim_write_spatial_maps:BadEventField", ...
            "event_field must be peak|pre|post.");
end

d = dir(outRoot);
movieDirs = d([d.isdir] & ~startsWith({d.name}, "."));
nOk = 0;
maps = {};
labels = {};

for i = 1:numel(movieDirs)
    movie = movieDirs(i).name;
    if strcmp(movie, "lab_meeting")
        continue
    end
    if ~isempty(onlyMovies) && ~any(contains(movie, onlyMovies))
        continue
    end
    stimIdx = find_stim_index(fnames, movie);
    if isempty(stimIdx)
        fprintf("SKIP  %s  (no stim entry)\n", movie);
        continue
    end
    extFile = find_ext_file(extRoot, movie);
    if isempty(extFile)
        fprintf("SKIP  %s  (no *_ext.mat)\n", movie);
        continue
    end

    fprintf("\n======== %s ========\n  ext: %s\n", movie, extFile);
    tw = h5read(extFile, "/temporal_weights");
    onset = round(double(anchors{stimIdx}(:)));
    onset = onset(isfinite(onset) & onset >= 1);
    if size(tw, 1) < size(tw, 2) && size(tw, 2) >= max(onset)
        tw = tw.';
    end
    T = double(tw);
    [~, nNeurons] = size(T);
    z = (T - mean(T, 1)) ./ max(std(T, 0, 1), 1e-8);
    clear T tw

    dsb = 1;
    try
        dsb = double(h5read(extFile, "/config/downsample_time_by"));
    catch
    end
    if ~isfinite(dsb) || dsb <= 0
        dsb = 1;
    end
    fs = fpsMovie / dsb;
    onsetTrace = max(1, round(onset / dsb));

    [amp, pval] = peri_amp_pval(z, fs, onsetTrace, periWin, baseWin, respWin);
    isExc = pval < alpha & amp > 0;
    isInh = pval < alpha & amp < 0;
    fprintf("  %d neurons | exc=%d (%.1f%%) inh=%d (%.1f%%)\n", ...
        nNeurons, sum(isExc), 100 * mean(isExc), sum(isInh), 100 * mean(isInh));

    [Sn, frameSize] = extract_load_ext_spatial(extFile, nNeurons, framePx);
    if frameSize(1) ~= frameSize(2)
        warning("extract_stim_write_spatial_maps:NonSquare", ...
            "FOV is %d x %d; using height=%d for reshape.", ...
            frameSize(1), frameSize(2), frameSize(1));
    end
    usePx = frameSize(1);

    sumimg = [];
    try
        sumimg = h5read(extFile, "/info/summary_image")';
    catch
    end

    outDir = fullfile(outRoot, movie);
    outFile = fullfile(outDir, "fig4_spatial_response.png");
    P = struct("Sn", Sn, "amp", amp, "isExc", isExc, "isInh", isInh, ...
        "framePx", usePx, "sumimg", sumimg, "outFile", outFile, "runId", movie);
    mapSmall = extract_stim_plot_spatial_fig(P);
    close(gcf);

    save(fullfile(outDir, "spatial_map.mat"), "amp", "pval", "isExc", "isInh", ...
        "mapSmall", "efLbl", "-v7.3");

    info = parse_rec_label(movie);
    maps{end + 1} = mapSmall; %#ok<AGROW>
    labels{end + 1} = info.label; %#ok<AGROW>
    copy_into_lab_meeting(outRoot, info.label, outFile);
    nOk = nOk + 1;
    clear z Sn amp pval
end

if nOk == 0
    error("extract_stim_write_spatial_maps:NoMovies", "No movies mapped under %s", outRoot);
end
write_cohort_gallery(outRoot, maps, labels, efLbl);
fprintf("\nWrote spatial maps for %d movies -> %s\n", nOk, outRoot);
end

function [amp, pval] = peri_amp_pval(z, fs, onset, periWin, baseWin, respWin)
pk = round(periWin(1) * fs):round(periWin(2) * fs);
tPeri = pk(:) / fs;
idxBase = tPeri >= baseWin(1) & tPeri < baseWin(2);
idxResp = tPeri >= respWin(1) & tPeri < respWin(2);
nFrames = size(z, 1);
nNeurons = size(z, 2);
valid = (onset + pk(1) >= 1) & (onset + pk(end) <= nFrames);
onsetV = onset(valid);
nT = numel(onsetV);
trialBase = zeros(nNeurons, nT);
trialResp = zeros(nNeurons, nT);
for k = 1:nT
    seg = double(z(onsetV(k) + pk, :));
    trialBase(:, k) = mean(seg(idxBase, :), 1).';
    trialResp(:, k) = mean(seg(idxResp, :), 1).';
end
amp = mean(trialResp - trialBase, 2);
pval = ones(nNeurons, 1);
for nrn = 1:nNeurons
    d = trialResp(nrn, :) - trialBase(nrn, :);
    if all(abs(d) < 1e-12)
        pval(nrn) = 1;
    else
        pval(nrn) = signrank(trialResp(nrn, :), trialBase(nrn, :));
    end
end
end

function idx = find_stim_index(fnames, movie)
idx = [];
for i = 1:numel(fnames)
    [~, stem, ~] = fileparts(fnames{i});
    if strcmp(stem, movie) || startsWith(movie, stem) || startsWith(stem, movie)
        idx = i;
        return
    end
end
end

function extFile = find_ext_file(extRoot, movie)
extFile = "";
d = dir(fullfile(extRoot, "*_ext.mat"));
if isempty(d)
    return
end
names = {d.name};
hit = find(startsWith(names, movie), 1);
if isempty(hit)
    hit = find(contains(names, movie), 1);
end
if ~isempty(hit)
    extFile = fullfile(extRoot, names{hit});
end
end

function copy_into_lab_meeting(outRoot, label, srcFile)
dst = fullfile(outRoot, "lab_meeting", "per_recording", label, "fig4_spatial_response.png");
dstDir = fileparts(dst);
if exist(dstDir, "dir") ~= 7
    mkdir(dstDir);
end
if exist(srcFile, "file") == 2
    copyfile(srcFile, dst);
end
end

function write_cohort_gallery(outRoot, maps, labels, efLbl)
avg = fullfile(outRoot, "lab_meeting", "cohort_average");
if exist(avg, "dir") ~= 7
    mkdir(avg);
end
n = numel(maps);
ncol = 4;
nrow = ceil(n / ncol);
cmapRb = [[linspace(0, 1, 128), ones(1, 128)]' ...
    [linspace(0, 1, 128), linspace(1, 0, 128)]' ...
    [ones(1, 128), linspace(1, 0, 128)]'];

figure("Position", [20, 20, 1680, 920], "Color", "w");
sgtitle({sprintf("Fig4  |  Poke response on FOV  (t=0 at %s)", efLbl), ...
    "Red = excited  |  blue = inhibited  |  each panel is one movie (not spatially aligned)"}, ...
    "FontSize", 12, "FontWeight", "bold");
for i = 1:n
    subplot(nrow, ncol, i);
    M = maps{i};
    climv = prctile(abs(M(M ~= 0)), 98);
    if isempty(climv) || ~isfinite(climv) || climv == 0
        climv = 1;
    end
    imagesc(M);
    colormap(gca, cmapRb);
    clim([-climv, climv]);
    axis image off;
    title(strrep(labels{i}, "_", "\_"), "FontSize", 8);
end
exportgraphics(gcf, fullfile(avg, "fig4_spatial_response.png"), "Resolution", 180);
close(gcf);
end

function info = parse_rec_label(movie)
info.mouse = "unknown";
info.session = "";
info.mirror = "";
info.power = "";
if contains(movie, "mirror1")
    info.mirror = "mirror1";
else
    info.mirror = "mirror2";
end
tok = regexp(movie, "s(\d+)", "tokens", "once");
if ~isempty(tok)
    info.session = sprintf('s%s', tok{1});
end
if contains(movie, "_m2c-") || contains(movie, "m201_")
    info.mouse = "m2c";
elseif contains(movie, "m3_")
    info.mouse = "m3_m1f";
elseif contains(movie, "_m1f-") || contains(movie, "_m1_")
    info.mouse = "m1f";
end
if strcmp(info.session, "s00")
    info.power = "8.38";
elseif strcmp(info.session, "s01")
    info.power = "20";
else
    info.power = "unknown";
end
pwr = info.power;
if strcmp(pwr, "unknown")
    pwr = "pwrNA";
else
    pwr = ["pwr" pwr];
end
info.label = sprintf("%s_%s_%s_%s", info.mouse, info.session, pwr, info.mirror);
end

function v = getp(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
