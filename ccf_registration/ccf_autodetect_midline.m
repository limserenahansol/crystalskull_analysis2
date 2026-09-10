function out = ccf_autodetect_midline(varargin)
%CCF_AUTODETECT_MIDLINE  Find the sagittal midline in every cached FOV, automatically.
%
% Each crystal-skull summary image shows both hemispheres as two blobs separated by a dark
% band - the superior sagittal sinus. That band is the midline: it fixes the rotation of the
% brain in the image and the ML = 0 line. Two of the four similarity degrees of freedom come
% free, without anyone clicking anything.
%
% What is left after this:
%   scale     known from the optics (um per pixel), not fitted - pass it in
%   AP offset one click along the midline, or one known landmark
%
%   out = ccf_autodetect_midline                 % all cached images + QC sheet
%   out = ccf_autodetect_midline("Only","poke")  % subset by label substring
%
% Writes cache/midline.csv and cache/midline_qc.png. Check the QC sheet before trusting any
% of it: the detector assumes exactly two hemisphere blobs and will misfire on a FOV where
% the window clips one hemisphere badly.

p = inputParser;
p.addParameter("Only", "", @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
only = string(p.Results.Only);

here = fileparts(mfilename("fullpath"));
cacheDir = fullfile(here, "cache");
files = dir(fullfile(cacheDir, "*.mat"));
files = files(~ismember(string({files.name}), ["fov_similarity.mat"]));
if only ~= ""
    files = files(contains(string({files.name}), only));
end
if isempty(files)
    error("ccf_autodetect_midline:NoCache", ...
        "No cached images in %s. Run ccf_cache_fov_images first.", cacheDir);
end

n = numel(files);
labels = strings(n, 1);
angleDeg = nan(n, 1);
midX = nan(n, 1);
midY = nan(n, 1);
sep = nan(n, 1);
quality = nan(n, 1);
imgs = cell(n, 1);
masks = cell(n, 1);

for k = 1:n
    [~, lab] = fileparts(files(k).name);
    labels(k) = string(lab);
    S = load(fullfile(cacheDir, files(k).name));
    img = double(S.img);
    if contains(lab, "mirror2")
        img = fliplr(img);      % same frame the registration will use
    end
    imgs{k} = img;
    try
        [angleDeg(k), midX(k), midY(k), sep(k), quality(k), masks{k}] = local_midline(img);
    catch ME
        fprintf(2, "  %-36s FAILED: %s\n", lab, ME.message);
    end
end

T = table(labels, angleDeg, midX, midY, sep, quality, ...
    VariableNames=["label", "angle_deg", "mid_x", "mid_y", "hemi_sep_px", "quality"]);
T = sortrows(T, "quality", "descend");
writetable(T, fullfile(cacheDir, "midline.csv"));

fprintf("\n%-36s %8s %8s %8s %9s %8s\n", "label", "angle", "mid_x", "mid_y", "sep_px", "quality");
for k = 1:height(T)
    fprintf("%-36s %8.2f %8.0f %8.0f %9.0f %8.2f\n", T.label(k), T.angle_deg(k), ...
        T.mid_x(k), T.mid_y(k), T.hemi_sep_px(k), T.quality(k));
end

atEdge = abs(abs(T.angle_deg) - 30) < 0.3;
good = T.quality > 0.12 & isfinite(T.angle_deg) & ~atEdge;
fprintf("\n%d/%d passed the quality bar\n", sum(good), height(T));
fprintf("angle across all: median %.2f deg, spread %.2f to %.2f\n", ...
    median(T.angle_deg(good), "omitnan"), min(T.angle_deg(good)), max(T.angle_deg(good)));

% ---- QC sheet
f = figure("Color", "w", "Position", [30 30 1750 980], "Visible", "off");
tl = tiledlayout(f, ceil(n / 6), 6, "TileSpacing", "compact", "Padding", "compact");
for k = 1:n
    ax = nexttile(tl);
    imshow(local_stretch(imgs{k}), "Parent", ax);
    hold(ax, "on");
    if isfinite(angleDeg(k))
        th = deg2rad(angleDeg(k));
        L = 1.2 * max(size(imgs{k}));
        d = [cos(th), sin(th)];
        p1 = [midX(k), midY(k)] - L * d;
        p2 = [midX(k), midY(k)] + L * d;
        plot(ax, [p1(1) p2(1)], [p1(2) p2(2)], "c-", "LineWidth", 1.6);
        plot(ax, midX(k), midY(k), "y+", "MarkerSize", 9, "LineWidth", 1.5);
    end
    col = "k";
    if quality(k) <= 0.12 || ~isfinite(quality(k)) || abs(abs(angleDeg(k)) - 30) < 0.3
        col = "r";
    end
    title(ax, sprintf("%s  q=%.2f  %.1f deg", labels(k), quality(k), angleDeg(k)), ...
        "FontSize", 6, "Interpreter", "none", "Color", col);
end
title(tl, "Auto-detected sagittal midline (cyan). Red title = low confidence, register by hand.");
exportgraphics(f, fullfile(cacheDir, "midline_qc.png"), "Resolution", 130);
close(f);
fprintf("wrote %s\n", fullfile(cacheDir, "midline_qc.png"));

out = T;
end


% ====================================================================== helpers
function [angDeg, cx, cy, sepPx, q, mask] = local_midline(img)
% The sagittal sinus is a dark line running the length of the brain, so it shows up as a
% trough in the row profile once the image is rotated to make it horizontal. Sweep the
% angle and take the rotation whose trough is deepest and narrowest.
%
% This replaces an earlier connected-component version that tried to split the two
% hemispheres and pick the perpendicular bisector of their centroids. That failed on 18/30
% FOVs because the threshold often split one hemisphere in two rather than separating them,
% giving nonsense angles like -45 or +88 degrees. Working from the projection uses the whole
% image and never has to make a binary decision about which blob is which.

SCALE = 512;
g0 = local_stretch(img);
scl = SCALE / max(size(g0));
g = imresize(g0, scl);
g = imgaussfilt(g, 2);

% brain mask, only so background does not dominate the row means
mask = imfill(imopen(g > graythresh(g) * 0.5, strel("disk", 3)), "holes");
if nnz(mask) < 0.05 * numel(mask)
    mask = g > prctile(g(:), 60);
end

angles = -30:0.25:30;
best = -inf;
angDeg = NaN;
bestRow = NaN;
bestProf = [];
for a = angles
    gr = imrotate(g, a, "bilinear", "crop");
    mr = imrotate(mask, a, "nearest", "crop");
    cnt = sum(mr, 2);
    prof = sum(gr .* mr, 2) ./ max(cnt, 1);
    prof(cnt < 0.25 * max(cnt)) = NaN;          % ignore rows barely inside the brain
    s = local_trough(prof);
    if s.score > best
        best = s.score;
        angDeg = a;
        bestRow = s.row;
        bestProf = prof;
    end
end

if ~isfinite(angDeg) || ~isfinite(bestRow)
    error("no midline trough found");
end

% trough row in the rotated frame -> a point on the midline in the original frame
c = (size(g) + 1) / 2;                           % imrotate 'crop' rotates about the centre
dy = bestRow - c(1);
th = deg2rad(angDeg);
% rotating the image by +a is rotating coordinates by -a; undo it
pRot = [0, dy];
pOrig = [pRot(1) * cosd(-angDeg) - pRot(2) * sind(-angDeg), ...
         pRot(1) * sind(-angDeg) + pRot(2) * cosd(-angDeg)];
cx = (c(2) + pOrig(1)) / scl;
cy = (c(1) + pOrig(2)) / scl;

% separation between the two hemisphere humps either side of the trough
sepPx = local_hump_sep(bestProf, bestRow) / scl;
q = best;
mask = imresize(mask, size(g0), "nearest");
end


function s = local_trough(prof)
% deepest narrow dip in the central half of the profile
s = struct("score", -inf, "row", NaN);
n = numel(prof);
valid = find(isfinite(prof));
if numel(valid) < 40
    return
end
lo = max(valid(1), round(n * 0.30));
hi = min(valid(end), round(n * 0.70));
if hi - lo < 20
    return
end
sm = movmean(prof, 5, "omitnan");
seg = sm(lo:hi);
[dip, i] = min(seg);
row = lo + i - 1;
% shoulders: the tissue on either side, a bit away from the dip
w = round(n * 0.06);
la = sm(max(valid(1), row - 3 * w):max(valid(1), row - w));
lb = sm(min(valid(end), row + w):min(valid(end), row + 3 * w));
shoulder = mean([mean(la, "omitnan"), mean(lb, "omitnan")], "omitnan");
if ~isfinite(shoulder) || shoulder <= 0 || ~isfinite(dip)
    return
end
% raw relative depth of the trough. An earlier version multiplied this by 3 and clipped to
% 1, which pinned all 30 FOVs at exactly 1.00 and made the score useless for ranking.
% Leave it unscaled so the values actually separate.
s.score = max((shoulder - dip) / shoulder, 0);
s.row = row;
end


function d = local_hump_sep(prof, row)
sm = movmean(prof, 5, "omitnan");
n = numel(sm);
a = sm(1:max(1, row - 1));
b = sm(min(n, row + 1):end);
[~, ia] = max(a);
[~, ib] = max(b);
if isempty(ia) || isempty(ib)
    d = NaN;
    return
end
d = (min(n, row + 1) + ib - 1) - ia;
end


function out = local_stretch(img)
v = img(isfinite(img));
if isempty(v)
    out = zeros(size(img));
    return
end
lo = prctile(v, 1);
hi = prctile(v, 99.5);
if hi <= lo
    hi = lo + 1;
end
out = min(max((img - lo) / (hi - lo), 0), 1);
out(~isfinite(out)) = 0;
end
