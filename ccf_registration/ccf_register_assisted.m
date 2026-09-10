function ccf_register_assisted(varargin)
%CCF_REGISTER_ASSISTED  Register every cached FOV with the fewest possible clicks.
%
% The midline is already found automatically (ccf_autodetect_midline), which fixes rotation
% and the ML = 0 line. That leaves the position along the midline and the scale. So:
%
%   with UmPerPx given   1 click per FOV  (the sinus confluence)
%   without it           2 clicks per FOV (anterior cortical pole, then the sinus
%                        confluence) and scale is fitted
%
% The skull is off in this prep so bregma and lambda are not visible. The clicks are on
% vascular landmarks instead - see local_midline_targets.
%
%   ccf_register_assisted                          % 2 clicks each, all 30
%   ccf_register_assisted("UmPerPx", 4.3)          % 1 click each
%   ccf_register_assisted("Only", "poke")          % just the face-poke FOVs
%   ccf_register_assisted("Redo", true)            % re-do ones already registered
%
% Resume-safe: anything already in transforms/ is skipped unless Redo is true.
% Writes transforms/<label>.mat and a QC png per FOV.

p = inputParser;
p.addParameter("Only", "", @(s) ischar(s) || isstring(s));
p.addParameter("UmPerPx", NaN, @isnumeric);
p.addParameter("Redo", false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
only = string(p.Results.Only);
umPerPx = p.Results.UmPerPx;

here = fileparts(mfilename("fullpath"));
addpath(here, fileparts(here));
cacheDir = fullfile(here, "cache");
outDir = fullfile(here, "transforms");
if exist(outDir, "dir") ~= 7
    mkdir(outDir);
end

midFile = fullfile(cacheDir, "midline.csv");
if exist(midFile, "file") ~= 2
    error("ccf_register_assisted:NoMidline", ...
        "Run ccf_cache_fov_images then ccf_autodetect_midline first.");
end
M = readtable(midFile, "Delimiter", ",", "TextType", "string");
if only ~= ""
    M = M(contains(M.label, only), :);
end
M = sortrows(M, "quality", "descend");

A = load(fullfile(here, "allen_dorsal_map.mat"));
resUm = double(A.res_um);
bregmaRC = double(A.bregma_rc) + 1;
lambdaRC = double(A.lambda_rc) + 1;

% The skull is off in this prep, so bregma and lambda sutures are NOT visible - an earlier
% version of this function asked for them, which cannot be done on these images. What IS
% clearly visible is the superior sagittal sinus and, at its posterior end, the Y where it
% splits into the two transverse sinuses (the confluence of sinuses). That confluence sits
% at approximately lambda. The anterior cortical pole on the midline is the other usable
% edge. Both coordinates are read from the atlas so they stay correct if the map is rebuilt.
[antPoleRC, confluenceRC] = local_midline_targets(A, bregmaRC, lambdaRC);

fprintf("\n%d FOV(s). ", height(M));
if isfinite(umPerPx)
    fprintf("Scale fixed at %.2f um/px, so 1 click each: the SINUS CONFLUENCE.\n", umPerPx);
else
    fprintf("%s", "Scale will be fitted, so 2 clicks each: the ANTERIOR cortical pole " + ...
        "on the" + newline + "midline, then the SINUS CONFLUENCE (the Y where the " + ...
        "sagittal sinus splits)." + newline);
end
fprintf("Close the window or press Escape at any point to stop; progress is kept.\n\n");

nDone = 0;
for k = 1:height(M)
    lab = M.label(k);
    dst = fullfile(outDir, lab + ".mat");
    if exist(dst, "file") == 2 && ~p.Results.Redo
        fprintf("  [%2d/%2d] %-36s already done, skipping\n", k, height(M), lab);
        continue
    end

    S = load(fullfile(cacheDir, lab + ".mat"));
    img = double(S.img);
    flipLR = contains(lab, "mirror2");
    if flipLR
        img = fliplr(img);
    end

    pts = local_click(img, M(k, :), umPerPx, lab, k, height(M));
    if isempty(pts)
        fprintf("\nstopped by user after %d newly registered FOV(s)\n", nDone);
        return
    end

    % atlas-side targets
    if isfinite(umPerPx)
        % one point: the sinus confluence. Rotation from the detected midline,
        % scale from the optics, position from the click.
        theta = deg2rad(M.angle_deg(k));
        s = umPerPx / resUm;                       % atlas px per FOV px
        R = s * [cos(theta), sin(theta); -sin(theta), cos(theta)];
        t = [confluenceRC(2); confluenceRC(1)] - R * pts(1, :)';
        tform = affinetform2d([R, t; 0 0 1]);
        fovPts = pts(1, :);
        ccfPts = [confluenceRC(2), confluenceRC(1)];
    else
        fovPts = pts;
        ccfPts = [antPoleRC(2), antPoleRC(1); confluenceRC(2), confluenceRC(1)];
        tform = fitgeotform2d(fovPts, ccfPts, "similarity");
    end

    sc = sqrt(abs(det(tform.A(1:2, 1:2))));
    fprintf("  [%2d/%2d] %-36s  %.2f um/px", k, height(M), lab, sc * resUm);
    if ~isfinite(umPerPx)
        fprintf("  (fitted)");
    end
    fprintf("\n");

    T = struct("tform", tform, "tform_type", "similarity", "fov_points", fovPts, ...
        "ccf_points", ccfPts, "res_um", resUm, "flipLR", flipLR, ...
        "fov_size", size(img), "um_per_px", sc * resUm, ...
        "midline_angle_deg", M.angle_deg(k), "midline_quality", M.quality(k), ...
        "map_file", string(fullfile(here, "allen_dorsal_map.mat")), ...
        "label", lab, "created", string(datetime("now")));
    save(dst, "-struct", "T");
    local_qc(img, A, tform, lab, fullfile(outDir, lab + "_qc.png"));

    % ccf_label_facepoke_n18 looks for n18/<recLabel>/ccf_transform.mat, so drop a copy
    % there for the face-poke recordings. Without this the labelling step silently skips
    % everything with "no transform".
    if startsWith(lab, "poke_")
        recLabel = extractAfter(lab, "poke_");
        recDir = fullfile(here, "n18", recLabel);
        if exist(recDir, "dir") ~= 7
            mkdir(recDir);
        end
        save(fullfile(recDir, "ccf_transform.mat"), "-struct", "T");
    end
    nDone = nDone + 1;
end

fprintf("\n%d newly registered, %d total in %s\n", nDone, ...
    numel(dir(fullfile(outDir, "*.mat"))), outDir);
fprintf("Check the _qc.png files, then run ccf_label_facepoke_n18 (or ccf_apply_to_cells).\n");
end


% ====================================================================== helpers
function pts = local_click(img, m, umPerPx, lab, k, n)
if isfinite(umPerPx)
    need = 1;
    what = "SINUS CONFLUENCE (Y junction, posterior end of the sagittal sinus)";
else
    need = 2;
    what = ["ANTERIOR cortical pole, on the midline", ...
            "SINUS CONFLUENCE (Y junction, posterior end)"];
end

f = figure("Name", sprintf("%d/%d  %s", k, n, lab), "NumberTitle", "off", ...
    "Color", "w", "Position", [60 60 950 950]);
imshow(local_stretch(img), []);
hold on;
% the detected midline, as a guide
th = deg2rad(m.angle_deg);
L = 1.3 * max(size(img));
d = [cos(th), sin(th)];
c = [m.mid_x, m.mid_y];
plot(c(1) + [-L L] * d(1), c(2) + [-L L] * d(2), "c-", "LineWidth", 1.2);
title({sprintf("%s   (midline %.1f deg, q %.2f)", lab, m.angle_deg, m.quality), ...
       "Click along the cyan midline. Escape to stop."}, "Interpreter", "none", "FontSize", 11);

pts = zeros(0, 2);
for i = 1:need
    xlabel(sprintf("click %s", what(i)), "FontSize", 13, "FontWeight", "bold");
    try
        [x, y] = ginput(1);
    catch
        pts = [];
        close(f);
        return
    end
    if isempty(x)
        pts = [];
        close(f);
        return
    end
    pts(end + 1, :) = [x, y]; %#ok<AGROW>
    plot(x, y, "y+", "MarkerSize", 16, "LineWidth", 2);
    text(x + 12, y, what(i), "Color", "y", "FontWeight", "bold");
end
close(f);
end


function local_qc(img, A, tform, lab, dst)
R = imref2d(size(img));
bw = imwarp(single(logical(A.boundary)), invert(tform), "nearest", "OutputView", R) > 0.5;
g = mat2gray(local_stretch(img));
rgb = repmat(g, 1, 1, 3);
c1 = rgb(:, :, 1); c2 = rgb(:, :, 2); c3 = rgb(:, :, 3);
c1(bw) = 1; c2(bw) = 0.2; c3(bw) = 0.2;
f = figure("Visible", "off", "Color", "w", "Position", [50 50 900 900]);
imshow(cat(3, c1, c2, c3));
title(lab, "Interpreter", "none");
exportgraphics(f, dst, "Resolution", 120);
close(f);
end


function [antPoleRC, confluenceRC] = local_midline_targets(A, bregmaRC, lambdaRC)
% Anterior cortical pole on the midline, read straight off the atlas, and the confluence
% of sinuses which we place at lambda.
lab = double(A.label_img);
mid = round(bregmaRC(2));
col = lab(:, mid);
acr = string(A.acronyms(:));
ids = double(A.ids(:));
olf = ["OLF", "MOB", "onl", "AOB"];
first = NaN;
for r = 1:numel(col)
    if col(r) == 0
        continue
    end
    j = find(ids == col(r), 1);
    if isempty(j) || any(acr(j) == olf)
        continue
    end
    first = r;
    break
end
if ~isfinite(first)
    first = bregmaRC(1) - 130;      % fallback, ~3.25 mm anterior
end
antPoleRC = [first, mid];
confluenceRC = [lambdaRC(1), mid];
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
