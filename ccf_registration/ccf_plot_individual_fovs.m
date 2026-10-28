function ccf_plot_individual_fovs(varargin)
%CCF_PLOT_INDIVIDUAL_FOVS  One Allen-superimposed figure per registered FOV.
%
% Each PNG shows the EXTRACT summary image (raw 1p FOV) next to that same
% image warped onto the Allen dorsal map, so vasculature and CCF borders
% can be compared. Face-poke recordings also get a poke-response panel.
%
%   ccf_plot_individual_fovs
%   ccf_plot_individual_fovs("Only", "poke")
%   ccf_plot_individual_fovs("Only", "m1f")
%
% Writes ccf_registration/figures/individual/<label>.png and two index
% montages in ccf_registration/figures/.

p = inputParser;
p.addParameter("Only", "", @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
only = string(p.Results.Only);

here = fileparts(mfilename("fullpath"));
addpath(here, fileparts(here));
outDir = fullfile(here, "figures", "individual");
if exist(outDir, "dir") ~= 7
    mkdir(outDir);
end

A = load(fullfile(here, "allen_dorsal_map.mat"));
Rfixed = imref2d(size(A.label_img));
atlasRgb = atlas_rgb(A);

d = dir(fullfile(here, "transforms", "*.mat"));
recs = string({d.name});
recs = extractBefore(recs, ".mat");
recs = recs(:);
if only ~= ""
    recs = recs(contains(recs, only));
end
if isempty(recs)
    error("ccf_plot_individual_fovs:None", ...
        "No transforms matched Only='%s' in %s.", only, fullfile(here, "transforms"));
end

fprintf("Writing %d individual FOV-on-Allen figures...\n", numel(recs));
thumbPoke = {};
namePoke = strings(0, 1);
thumbMel = {};
nameMel = strings(0, 1);

for i = 1:numel(recs)
    lab = recs(i);
    cacheFile = fullfile(here, "cache", lab + ".mat");
    if exist(cacheFile, "file") ~= 2
        fprintf("  skip %s (no cached FOV)\n", lab);
        continue
    end

    T = load(fullfile(here, "transforms", lab + ".mat"));
    C = load(cacheFile, "img");
    img = double(C.img);
    if T.flipLR
        img = fliplr(img);
    end
    fov = local_stretch(img);
    mask = brain_oval_mask(fov);

    fov512 = imresize(fov, [512, 512]);
    mask512 = imresize(double(mask), [512, 512]);
    t512 = scale_tform(T.tform, T.fov_size(1) / 512);
    brain = double(A.label_img > 0);
    warpedMask = imwarp(mask512, t512, "OutputView", Rfixed, "SmoothEdges", true) .* brain;
    warpedAnat = imwarp(fov512 .* mask512, t512, "OutputView", Rfixed, "SmoothEdges", true) .* brain;
    onAtlas = overlay_fov_on_atlas(atlasRgb, warpedAnat, warpedMask);
    onAtlas = paint_borders(onAtlas, A.boundary);

    nShow = 1024;
    fovShow = imresize(fov, [nShow, nShow]);
    tShow = scale_tform(T.tform, T.fov_size(1) / nShow);
    rawRgb = overlay_borders_on_fov(fovShow, A, tShow);
    umShow = T.um_per_px * (T.fov_size(1) / nShow);
    resp = [];
    if startsWith(lab, "poke_")
        resp = warp_poke_response(lab, T, Rfixed);
        if ~isempty(resp)
            resp(warpedMask < 0.35) = 0;
        end
    end

    outFile = fullfile(outDir, lab + ".png");
    write_one_figure(lab, T, rawRgb, onAtlas, resp, umShow, outFile);
    fprintf("  [%2d/%2d] %s\n", i, numel(recs), lab);

    thumb = imresize(onAtlas, [264, 228]);
    if startsWith(lab, "poke_")
        thumbPoke{end+1} = thumb; %#ok<AGROW>
        namePoke(end+1, 1) = lab; %#ok<AGROW>
    else
        thumbMel{end+1} = thumb; %#ok<AGROW>
        nameMel(end+1, 1) = lab; %#ok<AGROW>
    end
end

figDir = fullfile(here, "figures");
if ~isempty(thumbPoke)
    write_index(thumbPoke, namePoke, "Each face-poke FOV on Allen CCF (raw imaging warped)", ...
        fullfile(figDir, "fig_ccf_individual_index_poke.png"));
end
if ~isempty(thumbMel)
    write_index(thumbMel, nameMel, "Each melittin FOV on Allen CCF (raw imaging warped)", ...
        fullfile(figDir, "fig_ccf_individual_index_mel.png"));
end
fprintf("Wrote individual figures to %s\n", outDir);
end

function write_one_figure(lab, T, rawRgb, onAtlas, resp, umShow, outFile)
hasResp = ~isempty(resp) && any(isfinite(resp(:)) & resp(:) ~= 0);
if hasResp
    figW = 1680;
    nPanel = 3;
else
    figW = 1180;
    nPanel = 2;
end
fig = figure("Color", "w", "Position", [20, 20, figW, 720], "Visible", "off");

ax1 = subplot(1, nPanel, 1);
imshow(rawRgb, "Parent", ax1);
hold(ax1, "on");
draw_fov_guides(ax1, size(rawRgb, 1), umShow);
title(ax1, {"Raw FOV  (EXTRACT summary)", "P \leftarrow    \rightarrow A   ·   red = Allen borders"}, ...
    "FontSize", 11, "FontWeight", "normal");

ax2 = subplot(1, nPanel, 2);
imshow(onAtlas, "Parent", ax2);
hold(ax2, "on");
draw_area_labels(ax2);
title(ax2, "Same FOV warped onto Allen CCF", "FontSize", 11, "FontWeight", "normal");

if hasResp
    climv = prctile(abs(resp(resp ~= 0 & isfinite(resp))), 95);
    if ~(isfinite(climv) && climv > 0)
        climv = 1;
    end
    rgb = signed_on_atlas(onAtlas, resp, climv);
    ax3 = subplot(1, nPanel, 3);
    imshow(rgb, "Parent", ax3);
    hold(ax3, "on");
    draw_area_labels(ax3);
    title(ax3, {"Poke response on the same FOV", "red = excited, blue = inhibited"}, ...
        "FontSize", 11, "FontWeight", "normal");
end

diceTxt = "";
if isfield(T, "dice") && isfinite(T.dice)
    diceTxt = sprintf("   ·   Dice %.3f", T.dice);
end
sgtitle(fig, sprintf("%s     %.2f µm/px%s", pretty_title(lab), T.um_per_px, diceTxt), ...
    "FontSize", 14, "FontWeight", "bold", "Interpreter", "none");
exportgraphics(fig, outFile, "Resolution", 170);
close(fig);
end

function write_index(thumbs, names, ttl, outFile)
n = numel(thumbs);
ncol = 6;
nrow = ceil(n / ncol);
fig = figure("Color", "w", "Position", [20, 20, 240*ncol, 280*nrow], "Visible", "off");
for i = 1:n
    ax = subplot(nrow, ncol, i);
    imshow(thumbs{i}, "Parent", ax);
    title(ax, char(pretty_short(names(i))), "FontSize", 8, "Interpreter", "none");
end
sgtitle(ttl, "FontSize", 14, "FontWeight", "bold");
exportgraphics(fig, outFile, "Resolution", 150);
close(fig);
fprintf("  %s\n", outFile);
end

function rgb = overlay_fov_on_atlas(atlasRgb, warpedAnat, warpedMask)
rgb = atlasRgb;
m = warpedMask > 0.25;
if ~any(m(:))
    return
end
anat = warpedAnat;
v = anat(m);
lo = prctile(v, 2);
hi = prctile(v, 98);
if hi <= lo
    hi = lo + 1;
end
anat = min(max((anat - lo) / (hi - lo), 0), 1);
w = 0.84 * min(max(warpedMask, 0), 1);
for c = 1:3
    sl = rgb(:, :, c);
    sl(m) = (1 - w(m)) .* sl(m) + w(m) .* anat(m);
    rgb(:, :, c) = sl;
end
end

function rgb = overlay_borders_on_fov(fov, A, tform)
R = imref2d(size(fov));
bw = imwarp(single(logical(A.boundary)), invert(tform), "nearest", "OutputView", R) > 0.5;
bw = bwmorph(bw, "thin");
rgb = repmat(fov, 1, 1, 3);
ch1 = rgb(:, :, 1);
ch2 = rgb(:, :, 2);
ch3 = rgb(:, :, 3);
ch1(bw) = 0.92;
ch2(bw) = 0.18;
ch3(bw) = 0.18;
rgb = cat(3, ch1, ch2, ch3);
end

function rgb = paint_borders(rgb, boundary)
b = logical(boundary);
ch1 = rgb(:, :, 1);
ch2 = rgb(:, :, 2);
ch3 = rgb(:, :, 3);
ch1(b) = 0.80;
ch2(b) = 0.16;
ch3(b) = 0.16;
rgb = cat(3, ch1, ch2, ch3);
end

function draw_fov_guides(ax, nPx, umPerPx)
text(ax, 0.03*nPx, 0.50*nPx, "P", "Color", [1, 0.92, 0.35], "FontSize", 13, ...
    "FontWeight", "bold", "HorizontalAlignment", "left");
text(ax, 0.97*nPx, 0.50*nPx, "A", "Color", [1, 0.92, 0.35], "FontSize", 13, ...
    "FontWeight", "bold", "HorizontalAlignment", "right");
barUm = 1000;
barPx = barUm / umPerPx;
x0 = 0.08 * nPx;
y0 = 0.93 * nPx;
plot(ax, [x0, x0 + barPx], [y0, y0], "-", "Color", [1, 0.92, 0.35], "LineWidth", 2.4);
text(ax, x0 + barPx/2, y0 - 0.035*nPx, "1 mm", "Color", [1, 0.92, 0.35], ...
    "FontSize", 10, "HorizontalAlignment", "center", "FontWeight", "bold");
end

function t2 = scale_tform(tform, sc)
t2 = affinetform2d(tform.A * [sc, 0, 0; 0, sc, 0; 0, 0, 1]);
end

function resp = warp_poke_response(lab, T, Rfixed)
resp = [];
recLabel = extractAfter(lab, "poke_");
jobs = ccf_n18_jobs("Only", recLabel);
if height(jobs) == 0
    return
end
sp = fullfile(char(jobs.src_folder(1)), "spatial_map.mat");
if exist(sp, "file") ~= 2
    return
end
S = load(sp, "mapSmall");
M = imresize(double(S.mapSmall), T.fov_size);
if T.flipLR
    M = fliplr(M);
end
M512 = imresize(M, [512, 512]);
t512 = scale_tform(T.tform, T.fov_size(1) / 512);
resp = imwarp(M512, t512, "OutputView", Rfixed, "SmoothEdges", true);
resp = imgaussfilt(resp, 2.5);
end

function mask = brain_oval_mask(enh)
level = graythresh(enh);
bw = enh > max(level, 0.48);
bw = imclose(bw, strel("disk", 20));
bw = bwareafilt(bw, 1);
bw = imfill(bw, "holes");
mask = imerode(bw, strel("disk", 4));
end

function rgb = signed_on_atlas(atlasRgb, mu, climv)
cmap = red_blue_cmap();
z = mu / climv;
z = min(max(z, -1), 1);
idx = round((z + 1) / 2 * (size(cmap, 1) - 1)) + 1;
idx(~isfinite(mu) | mu == 0) = 1;
overlay = reshape(cmap(idx(:), :), [size(mu), 3]);
ok = isfinite(mu) & abs(mu) > 0.18 * climv;
w = 0.78 * ok;
rgb = atlasRgb;
for c = 1:3
    sl = atlasRgb(:, :, c);
    ov = overlay(:, :, c);
    sl(ok) = (1 - w(ok)) .* sl(ok) + w(ok) .* ov(ok);
    rgb(:, :, c) = sl;
end
end

function draw_area_labels(ax)
keep = ["MOs", "MOp", "SSp-bfd", "SSp-n", "SSp-m", "VISp", "RSPd", "AUDp"];
try
    P = ccf_area_centroid(keep, "mapFile", fullfile(fileparts(mfilename("fullpath")), ...
        "allen_dorsal_map.mat"));
    for i = 1:height(P)
        text(ax, P.x(i), P.y(i), char(P.acronym(i)), "FontSize", 7, ...
            "Color", [0.05, 0.05, 0.05], "HorizontalAlignment", "center", ...
            "FontWeight", "bold", "Interpreter", "none");
    end
catch
end
end

function bd = atlas_rgb(A)
d = double(A.surface_depth_um);
m = isfinite(d) & d > 0;
g = zeros(size(d));
if any(m(:))
    vv = d(m);
    lo = prctile(vv, 1);
    hi = prctile(vv, 99.5);
    g(m) = 1 - min(max((d(m) - lo) / max(hi - lo, eps), 0), 1);
end
bd = repmat(g, 1, 1, 3);
bd = paint_borders(bd, A.boundary);
end

function out = local_stretch(img)
v = img(isfinite(img));
if isempty(v)
    out = zeros(size(img));
    return
end
lo = prctile(v, 2);
hi = prctile(v, 97);
if hi <= lo
    hi = lo + 1;
end
out = min(max((img - lo) / (hi - lo), 0), 1);
out(~isfinite(out)) = 0;
out = out .^ 0.5;
end

function cmap = red_blue_cmap()
halfN = 128;
cmap = [[linspace(0, 1, halfN), ones(1, halfN)]', ...
    [linspace(0, 1, halfN), linspace(1, 0, halfN)]', ...
    [ones(1, halfN), linspace(1, 0, halfN)]'];
end

function s = pretty_title(lab)
lab = string(lab);
if startsWith(lab, "poke_")
    s = "Face poke   " + replace(extractAfter(lab, "poke_"), "_", "   ");
elseif startsWith(lab, "mel_")
    rest = replace(extractAfter(lab, "mel_"), "_", "   ");
    if contains(lab, "1014m0")
        s = "Melittin + morphine   " + rest;
    else
        s = "Melittin   " + rest;
    end
else
    s = lab;
end
end

function s = pretty_short(lab)
lab = string(lab);
if startsWith(lab, "poke_")
    s = extractAfter(lab, "poke_");
else
    s = extractAfter(lab, "mel_");
end
end
