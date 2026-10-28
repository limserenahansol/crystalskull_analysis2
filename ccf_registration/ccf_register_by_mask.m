function out = ccf_register_by_mask(varargin)
%CCF_REGISTER_BY_MASK  Register each FOV by matching its tissue silhouette to the atlas cortex.
%
% Why this replaces the click-based route
% ---------------------------------------
% Two clicked points carry four numbers. The tissue outline carries thousands, and it is
% already segmented. Clicking also failed here twice for a concrete reason: with the skull
% off there is no bregma, and the anterior cortical pole sits on the frontal curve where a
% top-down camera cannot see it, so the "anterior" click landed on the window edge and the
% fitted scale absorbed the error (5.2-12.8 um/px across recordings, and 25% of cells
% labelled cerebellum in a window that contains no cerebellum).
%
% What this does instead
%   rotation   from the auto-detected sagittal midline, plus the 90 deg between the FOV
%              frame (A-P horizontal) and the atlas frame (A-P vertical)
%   scale      + translation from matching the tissue mask to the atlas cortex mask,
%              refined by maximising Dice overlap
%   A-P sense  resolved by trying both and keeping the better overlap, since the cortex
%              silhouette is asymmetric front to back
%
%   out = ccf_register_by_mask                  % all cached FOVs
%   out = ccf_register_by_mask("Only","poke")
%   out = ccf_register_by_mask("DryRun",true)   % report, write nothing
%
% Writes transforms/<label>.mat (backing up any existing one to pre_mask/) and a QC png.

p = inputParser;
p.addParameter("Only", "", @(s) ischar(s) || isstring(s));
p.addParameter("DryRun", false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
only = string(p.Results.Only);
dryRun = logical(p.Results.DryRun);

here = fileparts(mfilename("fullpath"));
addpath(here);
cacheDir = fullfile(here, "cache");
tDir = fullfile(here, "transforms");
bDir = fullfile(tDir, "pre_mask");
if exist(tDir, "dir") ~= 7, mkdir(tDir); end

M = readtable(fullfile(cacheDir, "midline.csv"), "Delimiter", ",", "TextType", "string");
if only ~= ""
    M = M(contains(M.label, only), :);
end

A = load(fullfile(here, "allen_dorsal_map.mat"));
resUm = double(A.res_um);
ctxMask = local_cortex_mask(A);

n = height(M);
lab = strings(n, 1); umpp = nan(n, 1); dice = nan(n, 1); flipAP = false(n, 1);

for k = 1:n
    lab(k) = M.label(k);
    f = fullfile(cacheDir, lab(k) + ".mat");
    if exist(f, "file") ~= 2, continue; end
    S = load(f);
    img = double(S.img);
    if contains(lab(k), "mirror2"), img = fliplr(img); end

    fovMask = local_tissue_mask(img);
    [tform, d, fl] = local_fit(fovMask, ctxMask, M.angle_deg(k));
    dice(k) = d;
    flipAP(k) = fl;
    umpp(k) = sqrt(abs(det(tform.A(1:2, 1:2)))) * resUm;

    apTxt = "normal";
    if fl
        apTxt = "flipped";
    end
    fprintf("  %-36s dice %.3f   %.2f um/px   AP %s\n", lab(k), d, umpp(k), apTxt);
    if ~dryRun
        T = struct("tform", tform, "tform_type", "similarity", "res_um", resUm, ...
            "flipLR", contains(lab(k), "mirror2"), "fov_size", size(img), ...
            "um_per_px", umpp(k), "dice", d, "ap_flipped", fl, ...
            "midline_angle_deg", M.angle_deg(k), "method", "mask", ...
            "fov_points", [], "ccf_points", [], ...
            "map_file", string(fullfile(here, "allen_dorsal_map.mat")), ...
            "label", lab(k), "created", string(datetime("now")));
        old = fullfile(tDir, lab(k) + ".mat");
        if exist(old, "file") == 2
            if exist(bDir, "dir") ~= 7, mkdir(bDir); end
            if exist(fullfile(bDir, lab(k) + ".mat"), "file") ~= 2
                copyfile(old, fullfile(bDir, lab(k) + ".mat"));
            end
        end
        save(old, "-struct", "T");
        local_qc(img, A, tform, lab(k), d, fullfile(tDir, lab(k) + "_qc.png"));
        if startsWith(lab(k), "poke_")
            rd = fullfile(here, "n18", extractAfter(lab(k), "poke_"));
            if exist(rd, "dir") ~= 7, mkdir(rd); end
            save(fullfile(rd, "ccf_transform.mat"), "-struct", "T");
        end
    end
end

out = table(lab, umpp, dice, flipAP, VariableNames=["label", "um_per_px", "dice", "ap_flipped"]);
fprintf("\num/px: median %.2f, range %.2f-%.2f (spread %.2fx)\n", ...
    median(umpp, "omitnan"), min(umpp), max(umpp), max(umpp) / min(umpp));
fprintf("dice : median %.3f, worst %.3f\n", median(dice, "omitnan"), min(dice));
if dryRun, fprintf("\nDRY RUN - nothing written.\n"); end
end


% ====================================================================== helpers
function m = local_cortex_mask(A)
acr = string(A.acronyms(:));
ids = double(A.ids(:));
keep = ["MOs","MOp","SSp","SSs","VISp","VISa","VISam","VISal","VISpm","VISrl","VISl", ...
        "RSPd","RSPv","RSPagl","AUDp","ACAd","TEa","VISpor","SSp-bfd","SSp-m","SSp-n", ...
        "SSp-ul","SSp-ll","SSp-tr","SSp-un","AUDd","AUDv","AUDpo","ORBm","PL","ILA","FRP"];
sel = false(size(ids));
for i = 1:numel(ids)
    a = acr(i);
    base = extractBefore(a + "-", "-");
    sel(i) = any(a == keep) || any(base == keep);
end
m = ismember(double(A.label_img), ids(sel));
m = imfill(m, "holes");
end


function m = local_tissue_mask(img)
v = img(isfinite(img));
lo = prctile(v, 1); hi = prctile(v, 99.5);
g = min(max((img - lo) / max(hi - lo, 1), 0), 1);
g = imgaussfilt(g, 6);
m = g > prctile(g(:), 55);
m = imopen(m, strel("disk", 9));
m = imfill(m, "holes");
cc = bwconncomp(m);
if cc.NumObjects > 1
    a = cellfun(@numel, cc.PixelIdxList);
    m(:) = false;
    m(cc.PixelIdxList{find(a == max(a), 1)}) = true;
    % keep any other blob at least half the size of the biggest (the other hemisphere)
    for i = 1:numel(a)
        if a(i) > 0.4 * max(a), m(cc.PixelIdxList{i}) = true; end
    end
end
end


function [best, bestDice, bestFlip] = local_fit(fovMask, ctxMask, angleDeg)
% coarse: align centroid + extent after rotating the FOV into the atlas frame,
% then a small local search on scale and translation to maximise Dice.
[fr, fc] = find(fovMask);
fovC = [mean(fc), mean(fr)];
[ar, ac] = find(ctxMask);
ctxC = [mean(ac), mean(ar)];

best = []; bestDice = -inf; bestFlip = false;
for flipAP = [false true]
    % FOV A-P is horizontal, atlas A-P is vertical -> 90 deg, sign sets which end is anterior
    base = -angleDeg + 90;
    if flipAP, base = base + 180; end
    th = deg2rad(base);
    Rot = [cos(th), -sin(th); sin(th), cos(th)];

    % initial scale from the extent of the rotated FOV mask vs the cortex mask
    P = ([fc, fr] - fovC) * Rot';
    s0 = (max(ar) - min(ar)) / max(max(P(:, 2)) - min(P(:, 2)), eps);

    for s = s0 * (0.8:0.05:1.25)
        Rs = s * Rot;
        for dx = -30:15:30
            for dy = -30:15:30
                t = ctxC' + [dx; dy] - Rs * fovC';
                T = affinetform2d([Rs, t; 0 0 1]);
                d = local_dice(fovMask, ctxMask, T);
                if d > bestDice
                    bestDice = d; best = T; bestFlip = flipAP;
                end
            end
        end
    end
end

% refine
for iter = 1:3
    step = 8 / iter;
    improved = true;
    while improved
        improved = false;
        Acur = best.A;
        for ds = [-0.02 0 0.02] / iter
            for dx = [-step 0 step]
                for dy = [-step 0 step]
                    Anew = Acur;
                    Anew(1:2, 1:2) = Acur(1:2, 1:2) * (1 + ds);
                    Anew(1, 3) = Acur(1, 3) + dx;
                    Anew(2, 3) = Acur(2, 3) + dy;
                    T = affinetform2d(Anew);
                    d = local_dice(fovMask, ctxMask, T);
                    if d > bestDice + 1e-5
                        bestDice = d; best = T; improved = true;
                    end
                end
            end
        end
    end
end
end


function d = local_dice(fovMask, ctxMask, T)
R = imref2d(size(ctxMask));
w = imwarp(fovMask, T, "nearest", "OutputView", R);
i = nnz(w & ctxMask);
d = 2 * i / max(nnz(w) + nnz(ctxMask), 1);
end


function local_qc(img, A, tform, lab, d, dst)
R = imref2d(size(img));
bw = imwarp(single(logical(A.boundary)), invert(tform), "nearest", "OutputView", R) > 0.5;
v = img(isfinite(img));
lo = prctile(v, 1); hi = prctile(v, 99.5);
g = min(max((img - lo) / max(hi - lo, 1), 0), 1);
c1 = g; c2 = g; c3 = g;
c1(bw) = 1; c2(bw) = 0.2; c3(bw) = 0.2;
f = figure("Visible", "off", "Color", "w", "Position", [50 50 900 900]);
imshow(cat(3, c1, c2, c3));
title(sprintf("%s   dice %.3f", lab, d), "Interpreter", "none");
exportgraphics(f, dst, "Resolution", 110);
close(f);
end
