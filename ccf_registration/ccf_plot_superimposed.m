function ccf_plot_superimposed()
%CCF_PLOT_SUPERIMPOSED  All registered FOVs drawn on the Allen dorsal CCF.
%
% Uses ccf_registration/transforms/*.mat (from ccf_register_assisted).
% Writes PNGs to ccf_registration/figures/.

here = fileparts(mfilename("fullpath"));
addpath(here, fileparts(here));
outDir = fullfile(here, "figures");
if exist(outDir, "dir") ~= 7
    mkdir(outDir);
end

A = load(fullfile(here, "allen_dorsal_map.mat"));
Rfixed = imref2d(size(A.label_img));
atlasRgb = atlas_rgb(A);
groups = readtable(fullfile(here, "cache", "fov_groups.csv"), ...
    "Delimiter", ",", "TextType", "string");

d = dir(fullfile(here, "transforms", "*.mat"));
recs = strings(0, 1);
for i = 1:numel(d)
    recs(end+1, 1) = string(d(i).name(1:end-4)); %#ok<AGROW>
end
recs = unique(recs);
fprintf("Warping %d registered FOVs onto Allen CCF...\n", numel(recs));

W = struct("label", {}, "cohort", {}, "mouse", {}, "mask", {}, ...
    "anatomy", {}, "resp", {});
for i = 1:numel(recs)
    lab = recs(i);
    T = load(fullfile(here, "transforms", lab + ".mat"));
    cacheFile = fullfile(here, "cache", lab + ".mat");
    if exist(cacheFile, "file") ~= 2
        fprintf("  skip %s (no cached FOV)\n", lab);
        continue
    end
    C = load(cacheFile, "img");
    img = double(C.img);
    if T.flipLR
        img = fliplr(img);
    end
    img512 = imresize(local_stretch(img), [512, 512]);
    t512 = scale_tform(T.tform, T.fov_size(1) / 512);
    mask = brain_oval_mask(img512);
    brain = double(A.label_img > 0);
    warpedMask = imwarp(double(mask), t512, "OutputView", Rfixed, "SmoothEdges", true) .* brain;
    warpedAnat = imwarp(img512 .* double(mask), t512, "OutputView", Rfixed, "SmoothEdges", true) .* brain;

    g = groups(groups.label == lab, :);
    if height(g) == 0
        cohort = "unknown";
    else
        cohort = g.cohort(1);
    end

    W(end+1).label = lab; %#ok<AGROW>
    W(end).cohort = cohort;
    W(end).mouse = mouse_from_label(lab);
    W(end).mask = warpedMask;
    W(end).anatomy = warpedAnat;
    W(end).resp = [];
    if startsWith(lab, "poke_")
        W(end).resp = warp_poke_response(lab, T, Rfixed);
        if ~isempty(W(end).resp)
            W(end).resp(W(end).mask < 0.35) = 0;
        end
    end
    fprintf("  [%2d/%2d] %s\n", i, numel(recs), lab);
end

if isempty(W)
    error("ccf_plot_superimposed:Nothing", "No warped FOVs. Check transforms/ and cache/.");
end

poke = W(string({W.cohort}) == "face_poke");
mel = W(string({W.cohort}) ~= "face_poke");

write_outline_figure(atlasRgb, A, W, "All registered FOVs on Allen CCF", ...
    fullfile(outDir, "fig_ccf_all_outlines.png"), true);
write_outline_figure(atlasRgb, A, poke, "Face-poke FOVs (n=18) on Allen CCF", ...
    fullfile(outDir, "fig_ccf_poke_outlines.png"), true);
write_outline_figure(atlasRgb, A, mel, "Melittin FOVs on Allen CCF", ...
    fullfile(outDir, "fig_ccf_melittin_outlines.png"), true);

write_coverage_figure(atlasRgb, A, W, "FOV coverage (all recordings)", ...
    fullfile(outDir, "fig_ccf_all_coverage.png"));
write_coverage_figure(atlasRgb, A, poke, "FOV coverage (face poke)", ...
    fullfile(outDir, "fig_ccf_poke_coverage.png"));

write_mean_anatomy(atlasRgb, A, poke, "Face-poke anatomy mean on Allen CCF", ...
    fullfile(outDir, "fig_ccf_poke_anatomy.png"));
write_mean_response(atlasRgb, A, poke, "Face-poke response mean on Allen CCF", ...
    fullfile(outDir, "fig_ccf_poke_response.png"));
write_by_mouse(atlasRgb, A, poke, fullfile(outDir, "fig_ccf_poke_by_mouse.png"));
write_gallery(atlasRgb, poke, fullfile(outDir, "fig_ccf_poke_gallery.png"));
write_hero(atlasRgb, A, poke, fullfile(outDir, "fig_ccf_poke_superimposed.png"));

fprintf("Wrote figures to %s\n", outDir);
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


function write_outline_figure(atlasRgb, A, W, ttl, outFile, perMouse)
if nargin < 6
    perMouse = true;
end
fig = figure("Color", "w", "Position", [20, 20, 900, 980], "Visible", "off");
ax = axes(fig);
imshow(atlasRgb, "Parent", ax);
hold(ax, "on");
draw_area_labels(ax, A);
used = strings(0, 1);
hLeg = gobjects(0, 1);
if perMouse
    mice = unique(string({W.mouse}), "stable");
    for i = 1:numel(mice)
        sub = W(string({W.mouse}) == mice(i));
        U = false(size(sub(1).mask));
        for k = 1:numel(sub)
            U = U | (sub(k).mask > 0.35);
        end
        U = bwareafilt(imfill(U, "holes"), 1);
        [col, name] = color_for(sub(1).mouse, sub(1).cohort);
        B = bwboundaries(U, "noholes");
        if isempty(B)
            continue
        end
        [~, big] = max(cellfun(@(b) size(b, 1), B));
        hp = plot(ax, B{big}(:, 2), B{big}(:, 1), "-", "Color", col, "LineWidth", 2.0);
        used(end+1, 1) = name; %#ok<AGROW>
        hLeg(end+1, 1) = hp; %#ok<AGROW>
    end
else
    for i = 1:numel(W)
        [col, name] = color_for(W(i).mouse, W(i).cohort);
        U = bwareafilt(imfill(W(i).mask > 0.35, "holes"), 1);
        B = bwboundaries(U, "noholes");
        if isempty(B)
            continue
        end
        [~, big] = max(cellfun(@(b) size(b, 1), B));
        hp = plot(ax, B{big}(:, 2), B{big}(:, 1), "-", "Color", col, "LineWidth", 1.2);
        if ~any(used == name)
            used(end+1, 1) = name; %#ok<AGROW>
            hLeg(end+1, 1) = hp; %#ok<AGROW>
        end
    end
end
title(ax, ttl, "FontSize", 14, "FontWeight", "bold");
if ~isempty(hLeg)
    legend(hLeg, cellstr(used), "Location", "southoutside", "Orientation", "horizontal", ...
        "Box", "off", "Interpreter", "none");
end
exportgraphics(fig, outFile, "Resolution", 180);
close(fig);
fprintf("  %s\n", outFile);
end


function write_hero(atlasRgb, A, W, outFile)
cov = zeros(size(A.label_img));
acc = zeros(size(A.label_img));
nn = zeros(size(A.label_img));
for i = 1:numel(W)
    m = W(i).mask > 0.35;
    cov = cov + double(m);
    if ~isempty(W(i).resp)
        acc(m) = acc(m) + W(i).resp(m);
        nn(m) = nn(m) + 1;
    end
end
mu = nan(size(acc));
ok = nn > 0;
mu(ok) = acc(ok) ./ nn(ok);
climv = prctile(abs(mu(ok)), 95);
if ~(isfinite(climv) && climv > 0)
    climv = 1;
end
rgb = signed_on_atlas(atlasRgb, mu, climv);

fig = figure("Color", "w", "Position", [20, 20, 1500, 780], "Visible", "off");
ax1 = subplot(1, 2, 1);
imshow(atlasRgb, "Parent", ax1);
hold(ax1, "on");
h = imagesc(ax1, cov);
h.AlphaData = 0.55 * (cov > 0);
colormap(ax1, hot);
colorbar(ax1);
draw_area_labels(ax1, A);
title(ax1, sprintf("FOV coverage  (n=%d, max %d)", numel(W), max(cov(:))), ...
    "FontSize", 13, "FontWeight", "bold");

ax2 = subplot(1, 2, 2);
imshow(rgb, "Parent", ax2);
hold(ax2, "on");
mice = unique(string({W.mouse}), "stable");
for i = 1:numel(mice)
    sub = W(string({W.mouse}) == mice(i));
    U = false(size(sub(1).mask));
    for k = 1:numel(sub)
        U = U | (sub(k).mask > 0.35);
    end
    U = bwareafilt(imfill(U, "holes"), 1);
    B = bwboundaries(U, "noholes");
    if isempty(B)
        continue
    end
    [~, big] = max(cellfun(@(b) size(b, 1), B));
    [col, ~] = color_for(sub(1).mouse, sub(1).cohort);
    plot(ax2, B{big}(:, 2), B{big}(:, 1), "-", "Color", col, "LineWidth", 2.0);
end
draw_area_labels(ax2, A);
title(ax2, "Mean poke response + per-mouse FOV", "FontSize", 13, "FontWeight", "bold");
sgtitle("Face-poke recordings superimposed on Allen CCF", "FontSize", 16, "FontWeight", "bold");
exportgraphics(fig, outFile, "Resolution", 180);
close(fig);
fprintf("  %s\n", outFile);
end


function write_coverage_figure(atlasRgb, A, W, ttl, outFile)
cov = zeros(size(A.label_img));
for i = 1:numel(W)
    cov = cov + double(W(i).mask > 0.35);
end
fig = figure("Color", "w", "Position", [20, 20, 900, 980], "Visible", "off");
ax = axes(fig);
imshow(atlasRgb, "Parent", ax);
hold(ax, "on");
h = imagesc(ax, cov);
h.AlphaData = 0.55 * (cov > 0);
colormap(ax, hot);
cb = colorbar(ax);
cb.Label.String = "number of FOVs";
draw_area_labels(ax, A);
title(ax, sprintf("%s   (max %d)", ttl, max(cov(:))), "FontSize", 14, "FontWeight", "bold");
exportgraphics(fig, outFile, "Resolution", 180);
close(fig);
fprintf("  %s\n", outFile);
end


function write_mean_anatomy(atlasRgb, A, W, ttl, outFile)
acc = zeros(size(A.label_img));
n = zeros(size(A.label_img));
for i = 1:numel(W)
    m = W(i).mask > 0.35;
    acc(m) = acc(m) + W(i).anatomy(m);
    n(m) = n(m) + 1;
end
mu = zeros(size(acc));
ok = n > 0;
mu(ok) = acc(ok) ./ n(ok);
fig = figure("Color", "w", "Position", [20, 20, 900, 980], "Visible", "off");
ax = axes(fig);
imshow(imfuse(atlasRgb, mu, "falsecolor"), "Parent", ax);
hold(ax, "on");
draw_area_labels(ax, A);
title(ax, ttl, "FontSize", 14, "FontWeight", "bold");
exportgraphics(fig, outFile, "Resolution", 180);
close(fig);
fprintf("  %s\n", outFile);
end


function write_mean_response(atlasRgb, A, W, ttl, outFile)
acc = zeros(size(A.label_img));
n = zeros(size(A.label_img));
for i = 1:numel(W)
    if isempty(W(i).resp)
        continue
    end
    m = W(i).mask > 0.35;
    acc(m) = acc(m) + W(i).resp(m);
    n(m) = n(m) + 1;
end
if ~any(n(:) > 0)
    fprintf("  no poke response maps, skip %s\n", outFile);
    return
end
mu = nan(size(acc));
ok = n > 0;
mu(ok) = acc(ok) ./ n(ok);
climv = prctile(abs(mu(ok)), 95);
if ~(isfinite(climv) && climv > 0)
    climv = 1;
end
rgb = signed_on_atlas(atlasRgb, mu, climv);
fig = figure("Color", "w", "Position", [20, 20, 900, 980], "Visible", "off");
ax = axes(fig);
imshow(rgb, "Parent", ax);
hold(ax, "on");
draw_area_labels(ax, A);
title(ax, {ttl, "red = excited, blue = inhibited (mean across recordings)"}, ...
    "FontSize", 13, "FontWeight", "bold");
exportgraphics(fig, outFile, "Resolution", 180);
close(fig);
fprintf("  %s\n", outFile);
end


function write_by_mouse(atlasRgb, A, W, outFile)
mice = unique(string({W.mouse}), "stable");
n = numel(mice);
fig = figure("Color", "w", "Position", [20, 20, 380*n, 860], "Visible", "off");
for i = 1:n
    sub = W(string({W.mouse}) == mice(i));
    ax = subplot(1, n, i);
    acc = zeros(size(A.label_img));
    nn = zeros(size(A.label_img));
    hold(ax, "on");
    imshow(atlasRgb, "Parent", ax);
    for k = 1:numel(sub)
        if ~isempty(sub(k).resp)
            m = sub(k).mask > 0.35;
            acc(m) = acc(m) + sub(k).resp(m);
            nn(m) = nn(m) + 1;
        end
        B = bwboundaries(sub(k).mask > 0.35, "noholes");
        [col, ~] = color_for(sub(k).mouse, sub(k).cohort);
        for b = 1:numel(B)
            plot(ax, B{b}(:, 2), B{b}(:, 1), "-", "Color", col, "LineWidth", 1.2);
        end
    end
    if any(nn(:) > 0)
        mu = nan(size(acc));
        ok = nn > 0;
        mu(ok) = acc(ok) ./ nn(ok);
        climv = prctile(abs(mu(ok)), 95);
        if ~(isfinite(climv) && climv > 0)
            climv = 1;
        end
        rgb = signed_on_atlas(atlasRgb, mu, climv);
        imshow(rgb, "Parent", ax);
        hold(ax, "on");
        for k = 1:numel(sub)
            B = bwboundaries(sub(k).mask > 0.35, "noholes");
            for b = 1:numel(B)
                plot(ax, B{b}(:, 2), B{b}(:, 1), "-", "Color", [0.1, 0.1, 0.1], "LineWidth", 0.8);
            end
        end
    end
    draw_area_labels(ax, A);
    title(ax, sprintf("%s   (%d FOVs)", mice(i), numel(sub)), ...
        "FontSize", 13, "Interpreter", "none");
end
sgtitle("Face-poke FOVs + mean poke response on Allen CCF", "FontSize", 15, "FontWeight", "bold");
exportgraphics(fig, outFile, "Resolution", 170);
close(fig);
fprintf("  %s\n", outFile);
end


function write_gallery(atlasRgb, W, outFile)
n = numel(W);
ncol = 6;
nrow = ceil(n / ncol);
fig = figure("Color", "w", "Position", [20, 20, 280*ncol, 260*nrow], "Visible", "off");
for i = 1:n
    ax = subplot(nrow, ncol, i);
    imshow(imfuse(atlasRgb, W(i).anatomy, "falsecolor"), "Parent", ax);
    title(ax, char(W(i).label), "FontSize", 7, "Interpreter", "none");
end
sgtitle("Each face-poke FOV warped onto Allen CCF (magenta = atlas, green = FOV)", ...
    "FontSize", 14);
exportgraphics(fig, outFile, "Resolution", 140);
close(fig);
fprintf("  %s\n", outFile);
end


function rgb = signed_on_atlas(atlasRgb, mu, climv)
cmap = red_blue_cmap();
z = mu / climv;
z = min(max(z, -1), 1);
idx = round((z + 1) / 2 * (size(cmap, 1) - 1)) + 1;
idx(~isfinite(mu)) = 1;
overlay = reshape(cmap(idx(:), :), [size(mu), 3]);
w = 0.72 * isfinite(mu);
rgb = atlasRgb;
for c = 1:3
    sl = atlasRgb(:, :, c);
    ov = overlay(:, :, c);
    sl(isfinite(mu)) = (1 - w(isfinite(mu))) .* sl(isfinite(mu)) + w(isfinite(mu)) .* ov(isfinite(mu));
    rgb(:, :, c) = sl;
end
b = edge(isfinite(mu), "sobel");
for c = 1:3
    sl = rgb(:, :, c);
    sl(b) = 0.15;
    rgb(:, :, c) = sl;
end
end


function draw_area_labels(ax, ~)
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


function [col, name] = color_for(mouse, cohort)
mouse = string(mouse);
cohort = string(cohort);
name = mouse;
if cohort == "face_poke"
    switch mouse
        case "m1f"
            col = [0.12, 0.40, 0.82];
        case "m2c"
            col = [0.10, 0.60, 0.32];
        case "m3"
            col = [0.88, 0.45, 0.08];
        otherwise
            col = [0.25, 0.25, 0.25];
    end
elseif cohort == "melittin_morphine"
    col = [0.78, 0.16, 0.38];
    name = mouse + " +morphine";
else
    col = [0.52, 0.18, 0.68];
    name = mouse + " melittin";
end
end


function m = mouse_from_label(lab)
lab = string(lab);
if startsWith(lab, "poke_m3")
    m = "m3";
elseif startsWith(lab, "poke_m1f")
    m = "m1f";
elseif startsWith(lab, "poke_m2c")
    m = "m2c";
elseif contains(lab, "1014m01b")
    m = "1014m01b";
elseif contains(lab, "1014m02w")
    m = "1014m02w";
elseif contains(lab, "CS1014-1b")
    m = "CS1014-1b";
elseif contains(lab, "CS1014-2w")
    m = "CS1014-2w";
elseif contains(lab, "CS0204-1b")
    m = "CS0204-1b";
elseif contains(lab, "CS0204-2w")
    m = "CS0204-2w";
else
    m = lab;
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
b = logical(A.boundary);
ch1 = bd(:, :, 1); ch2 = bd(:, :, 2); ch3 = bd(:, :, 3);
ch1(b) = 0.75; ch2(b) = 0.20; ch3(b) = 0.20;
bd = cat(3, ch1, ch2, ch3);
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
