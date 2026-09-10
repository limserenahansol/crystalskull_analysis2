function T = ccf_register_fov(fovImage, varargin)
%CCF_REGISTER_FOV  Register a top-view crystal-skull FOV to the Allen CCF dorsal surface.
%
% The atlas side comes from allen_dorsal_map.mat, built by make_allen_dorsal_map.py from the
% BrainGlobe allen_mouse_25um volume: a top-down projection of the first labelled voxel in each
% (AP, ML) column, with cortical layers collapsed to their parent area. 25 um per atlas pixel.
%
%   T = ccf_register_fov(img)                       % bregma/lambda mode (2 clicks), the default
%   T = ccf_register_fov(img, "mode", "landmarks")  % cpselect, any number of matched points
%   T = ccf_register_fov(img, "mode", "points", "fovPoints", P, "ccfPoints", Q)   % scripted
%
% fovImage can be
%   - a numeric matrix (e.g. the EXTRACT summary image)
%   - a path to an EXTRACT *_ds_ext.h5 (reads /info/summary_image)
%   - a path to any image file
%
% Name/value options
%   mode        "bregma" (default) | "landmarks" | "points"
%   transform   "similarity" (default) | "affine" | "projective"
%   flipLR      false (default). Set true for mirror2, matching the convention in
%               extract_melittin_lab_meeting.m and the face-poke scripts.
%   orient      "Pleft_Aright" (default, matches lab maps) or "A_up" (native CCF)
%   mapFile     path to allen_dorsal_map.mat
%   outFile     where to save the transform (default ccf_transform_<timestamp>.mat)
%   fovPoints   Nx2 [x y] in FOV pixels          (mode "points")
%   ccfPoints   Nx2 [x y] in atlas pixels        (mode "points")
%
% Returns a struct with the fitted transform plus everything needed to label cells later.
% Feed it to ccf_apply_to_cells.m.
%
% Which mode to use
%   "bregma"     fastest and usually good enough: two clicks, similarity transform. Accuracy is
%                limited by how well you can see bregma/lambda, typically 200-400 um.
%   "landmarks"  better: add the sagittal sinus along the midline and the two transverse sinuses,
%                which are the clearest features in a top-view 1p image. 4-6 points, affine.
%   "points"     for functional registration - pass the peak coordinates of sensory maps
%                (visual flash -> VISp, whisker -> SSp-bfd, forelimb -> SSp-ul) together with
%                their atlas coordinates from ccf_area_centroid.m. This is the most accurate
%                option and the one to use if you ever run a mapping session.

arguments
    fovImage
end
arguments (Repeating)
    varargin
end

p = inputParser;
p.addParameter("mode", "bregma", @(s) ismember(string(s), ["bregma", "landmarks", "points"]));
p.addParameter("transform", "similarity", @(s) ismember(string(s), ["similarity", "affine", "projective"]));
p.addParameter("flipLR", false, @(x) islogical(x) || isnumeric(x));
p.addParameter("orient", "Pleft_Aright", @(s) ismember(string(s), ["Pleft_Aright", "A_up"]));
p.addParameter("mapFile", fullfile(fileparts(mfilename("fullpath")), "allen_dorsal_map.mat"), @(s) ischar(s) || isstring(s));
p.addParameter("outFile", "", @(s) ischar(s) || isstring(s));
p.addParameter("fovPoints", [], @isnumeric);
p.addParameter("ccfPoints", [], @isnumeric);
p.addParameter("saveQc", true, @(x) islogical(x) || isnumeric(x));
p.addParameter("closeQc", false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;
mode = string(opt.mode);
orient = string(opt.orient);

% ---------------------------------------------------------------- atlas
if exist(opt.mapFile, "file") ~= 2
    error("ccf_register_fov:NoMap", ...
        "Atlas map not found: %s\nRun make_allen_dorsal_map.py first.", opt.mapFile);
end
A = load(opt.mapFile);
A.label_img = double(A.label_img);
resUm = double(A.res_um);
% saved 0-based from Python
bregmaRC = double(A.bregma_rc) + 1;
lambdaRC = double(A.lambda_rc) + 1;

% ---------------------------------------------------------------- FOV image
fov = local_load_image(fovImage);
if opt.flipLR
    fov = fliplr(fov);
end
fovDisp = local_stretch(fov);

atlasNative = local_atlas_backdrop(A, bregmaRC, lambdaRC);
useRot = orient == "Pleft_Aright";
nAP = size(atlasNative, 1);
if useRot
    atlasDisp = imrotate(atlasNative, -90);   % P left, A right, left hemisphere on top
else
    atlasDisp = atlasNative;
end

% ---------------------------------------------------------------- get matched points
switch mode
    case "bregma"
        [fovPts, ccfPts] = local_pick_bregma(fovDisp, bregmaRC, lambdaRC);
        if size(fovPts, 1) < 2
            error("ccf_register_fov:Cancelled", "Need bregma and lambda; got %d point(s).", size(fovPts, 1));
        end
        tformType = "similarity";   % 2 points can only support similarity

    case "landmarks"
        if exist("cpselect", "file") ~= 2
            error("ccf_register_fov:NoIPT", ...
                "cpselect needs the Image Processing Toolbox. Use mode 'bregma' or 'points'.");
        end
        if useRot
            fprintf(['\nAtlas is rotated to match the lab maps: P left, A right, left hemisphere up.\n' ...
                     'Click the SAME anatomy on FOV (left) and atlas (right), then File > Close.\n' ...
                     'Best points: 2-3 along the midline (sagittal sinus), plus bregma if you can\n' ...
                     'see it, and the posterior sinus junction. Skip olfactory bulbs unless they\n' ...
                     'are actually in the crystal-skull window.\n\n']);
        else
            fprintf(['\nPick matching points in both images, then File > Close to continue.\n' ...
                     'Atlas here is A-up (standard CCF). Good landmarks: midline, bregma, lambda.\n\n']);
        end
        [fovPts, ccfClicked] = cpselect(fovDisp, atlasDisp, "Wait", true);
        if useRot
            ccfPts = local_rot_xy_to_native(ccfClicked, nAP);
        else
            ccfPts = ccfClicked;
        end
        tformType = string(opt.transform);

    case "points"
        fovPts = opt.fovPoints;
        ccfPts = opt.ccfPoints;
        if isempty(fovPts) || size(fovPts, 1) ~= size(ccfPts, 1)
            error("ccf_register_fov:BadPoints", "fovPoints and ccfPoints must be Nx2 and the same height.");
        end
        tformType = string(opt.transform);

    otherwise
        error("ccf_register_fov:BadMode", "mode must be bregma, landmarks, or points.");
end

nPts = size(fovPts, 1);
minPts = struct("similarity", 2, "affine", 3, "projective", 4);
if nPts < minPts.(tformType)
    error("ccf_register_fov:TooFewPoints", ...
        "%s needs at least %d points; you gave %d.", tformType, minPts.(tformType), nPts);
end

% ---------------------------------------------------------------- fit
tform = fitgeotform2d(fovPts, ccfPts, tformType);

% residuals, in microns on the atlas
pred = transformPointsForward(tform, fovPts);
residPx = sqrt(sum((pred - ccfPts).^2, 2));
residUm = residPx * resUm;

fprintf("\n%s fit on %d point(s)\n", tformType, nPts);
fprintf("  residual per point (um): %s\n", strjoin(compose("%.0f", residUm), ", "));
fprintf("  RMS residual: %.0f um\n", sqrt(mean(residUm.^2)));
if tformType == "similarity"
    sc = sqrt(abs(det(tform.A(1:2, 1:2))));
    fprintf("  implied scale: %.3f atlas px per FOV px  =  %.2f um per FOV px\n", sc, sc * resUm);
end
if nPts <= minPts.(tformType)
    fprintf("  NOTE: exactly determined, so the residual is 0 by construction and means nothing.\n");
    fprintf("        Add points if you want a real error estimate.\n");
end

% ---------------------------------------------------------------- pack
T = struct();
T.tform = tform;               % FOV pixels -> native atlas pixels (A-up)
T.tform_type = tformType;
T.fov_points = fovPts;
T.ccf_points = ccfPts;
T.resid_um = residUm;
T.rms_um = sqrt(mean(residUm.^2));
T.res_um = resUm;
T.flipLR = logical(opt.flipLR);
T.orient = orient;
T.fov_size = size(fov);
T.map_file = string(opt.mapFile);
T.bregma_rc = bregmaRC;
T.created = string(datetime("now"));

outFile = string(opt.outFile);
if outFile == ""
    outFile = fullfile(fileparts(mfilename("fullpath")), ...
        "ccf_transform_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".mat");
end
save(outFile, "-struct", "T");
fprintf("  saved -> %s\n", outFile);
T.out_file = outFile;

% ---------------------------------------------------------------- QC overlay
qcPng = "";
if logical(opt.saveQc)
    qcPng = char(replace(outFile, ".mat", "_qc.png"));
end
local_qc_figure(fovDisp, A, tform, fovPts, ccfPts, resUm, qcPng, logical(opt.closeQc));
T.qc_file = string(qcPng);
end


% ====================================================================== helpers
function img = local_load_image(src)
if isnumeric(src)
    img = double(src);
    return
end
src = char(src);
[~, ~, ext] = fileparts(src);
switch lower(ext)
    case {".h5", ".hdf5"}
        img = double(h5read(src, "/info/summary_image")');
    case ".mat"
        S = load(src);
        f = fieldnames(S);
        img = double(S.(f{1}));
    otherwise
        img = double(imread(src));
        if ndims(img) == 3
            img = mean(img, 3);
        end
end
end


function out = local_stretch(img)
% contrast stretch to [0 1]; the crystal-skull summary images are very low contrast
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
out = (img - lo) / (hi - lo);
out(~isfinite(out)) = 0;
out = min(max(out, 0), 1);
end


function bd = local_atlas_backdrop(A, bregmaRC, lambdaRC)
% grey brain with red-ish area borders, as a single RGB image for cpselect
d = double(A.surface_depth_um);
m = isfinite(d) & d > 0;
g = zeros(size(d));
if any(m(:))
    g(m) = 1 - local_stretch(d(m) * 0 + d(m));
end
bd = repmat(g, 1, 1, 3);
b = logical(A.boundary);
ch1 = bd(:, :, 1); ch2 = bd(:, :, 2); ch3 = bd(:, :, 3);
ch1(b) = 1; ch2(b) = 0.2; ch3(b) = 0.2;
bd = cat(3, ch1, ch2, ch3);
if nargin >= 3
    bd = local_draw_cross(bd, bregmaRC(2), bregmaRC(1), [0 1 1]);
    bd = local_draw_cross(bd, lambdaRC(2), lambdaRC(1), [1 0.85 0.1]);
end
end


function img = local_draw_cross(img, x, y, col)
[h, w, ~] = size(img);
x = round(x); y = round(y);
hw = 7;
rr = max(1, y-1):min(h, y+1);
cc = max(1, x-hw):min(w, x+hw);
for k = 1:3
    sl = img(:, :, k);
    sl(rr, cc) = col(k);
    img(:, :, k) = sl;
end
rr = max(1, y-hw):min(h, y+hw);
cc = max(1, x-1):min(w, x+1);
for k = 1:3
    sl = img(:, :, k);
    sl(rr, cc) = col(k);
    img(:, :, k) = sl;
end
end


function nativePts = local_rot_xy_to_native(rotPts, nAP)
% Clicks on imrotate(atlas, -90): x = new col, y = new row.
% Inverse of (old_row, old_col) -> (old_col, nAP+1-old_row).
x = rotPts(:, 1);
y = rotPts(:, 2);
nativeCol = y;
nativeRow = nAP + 1 - x;
nativePts = [nativeCol, nativeRow];
end


function [fovPts, ccfPts] = local_pick_bregma(fovDisp, bregmaRC, lambdaRC)
f = figure("Name", "Click bregma, then lambda", "NumberTitle", "off", "Color", "w");
imshow(fovDisp, []);
colormap(gray);
title({"Click BREGMA, then LAMBDA.", ...
       "Zoom with the toolbar first if you need to; press Enter when both are placed."}, ...
      "FontSize", 12);
hold on;
pts = zeros(0, 2);
for k = 1:2
    lbl = "bregma";
    if k == 2
        lbl = "lambda";
    end
    xlabel(sprintf("click %s", lbl), "FontSize", 12, "FontWeight", "bold");
    [x, y] = ginput(1);
    if isempty(x)
        break
    end
    pts(end + 1, :) = [x, y]; %#ok<AGROW>
    plot(x, y, "c+", "MarkerSize", 14, "LineWidth", 2);
    text(x + 8, y, lbl, "Color", "c", "FontWeight", "bold");
end
close(f);
fovPts = pts;
% atlas side: [x y] = [col row]
ccfPts = [bregmaRC(2), bregmaRC(1); lambdaRC(2), lambdaRC(1)];
ccfPts = ccfPts(1:size(fovPts, 1), :);
end


function local_qc_figure(fovDisp, A, tform, fovPts, ccfPts, resUm, qcPng, doClose)
% warp the atlas boundaries back into FOV space and lay them on the image
inv = invert(tform);
R = imref2d(size(fovDisp));
bw = imwarp(single(logical(A.boundary)), inv, "nearest", "OutputView", R) > 0.5;

f = figure("Name", "CCF registration QC", "Color", "w", "Position", [80 80 1400 650]);

subplot(1, 2, 1);
rgb = repmat(mat2gray(fovDisp), 1, 1, 3);
c1 = rgb(:, :, 1); c2 = rgb(:, :, 2); c3 = rgb(:, :, 3);
c1(bw) = 1; c2(bw) = 0.15; c3(bw) = 0.15;
imshow(cat(3, c1, c2, c3));
hold on;
plot(fovPts(:, 1), fovPts(:, 2), "c+", "MarkerSize", 14, "LineWidth", 2);
title("Atlas borders warped onto your FOV", "FontSize", 12);

subplot(1, 2, 2);
imshow(local_atlas_backdrop(A));
hold on;
pred = transformPointsForward(tform, fovPts);
plot(ccfPts(:, 1), ccfPts(:, 2), "co", "MarkerSize", 11, "LineWidth", 2);
plot(pred(:, 1), pred(:, 2), "yx", "MarkerSize", 11, "LineWidth", 2);
% FOV outline pushed into atlas space
[h, w] = size(fovDisp);
corners = transformPointsForward(tform, [1 1; w 1; w h; 1 h; 1 1]);
plot(corners(:, 1), corners(:, 2), "g-", "LineWidth", 2);
legend({"clicked", "predicted", "FOV outline"}, "TextColor", "k", "Location", "southoutside", ...
    "Orientation", "horizontal");
title(sprintf("Atlas space (A-up)  -  RMS residual %.0f um", ...
    sqrt(mean((sqrt(sum((pred - ccfPts).^2, 2)) * resUm).^2))), ...
    "FontSize", 12);
drawnow;
if nargin >= 7 && strlength(string(qcPng)) > 0
    exportgraphics(f, qcPng, "Resolution", 140);
    fprintf("  QC -> %s\n", qcPng);
end
if nargin >= 8 && doClose
    close(f);
end
end
