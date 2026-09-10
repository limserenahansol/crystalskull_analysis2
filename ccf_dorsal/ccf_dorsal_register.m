function S = ccf_dorsal_register(fov, varargin)
%CCF_DORSAL_REGISTER  Manual landmark warp of a top-view FOV to dorsal CCF.
%
%   S = ccf_dorsal_register(fov)
%   S = ccf_dorsal_register(fov, 'FlipLR', true, 'OutDir', outFolder)
%
% fov can be:
%   - a 2-D image
%   - a PNG/TIF path
%   - an EXTRACT *_ext.mat (uses /info/summary_image)
%   - a session folder (resolved with extract_resolve_extract_files)
%
% This is the right geometry for crystal-skull / 1p large-FOV imaging:
% a dorsal (top-down) Allen map, not a 3-D CCF volume slice.
%
% Better (recommended) two-step:
%   1) 1p skull-off photo  ->  Allen dorsal map  (bregma, lambda, midline, OB)
%   2) mesoscope FOV       ->  same 1p photo     (shared vessels)
%   Then compose:  ccf_dorsal_compose(S_fov_to_1p, S_1p_to_ccf)
%
% Suggested clicks (same order on FOV and atlas):
%   1 left olfactory bulb   2 right olfactory bulb
%   3 bregma (midline)      4 lambda / caudal RSC
%   then matching vessels or hemisphere edges (6–10 points total).
%
% Name-value:
%   FlipLR     - flip left-right (use for mirror2, P-left/A-right)
%   Photo1p    - optional 1p photo; if set, register FOV->1p first
%   AtlasPng   - optional Allen top-down PNG (overrides built-in sketch)
%   Transform  - 'similarity' | 'affine' (default) | 'projective'
%   Orient     - 'Pleft_Aright' (default, matches lab maps) or 'A_up'
%   OutDir     - folder for QC PNG + registration.mat

p = inputParser;
addParameter(p, "FlipLR", false);
addParameter(p, "Photo1p", "", @(s) ischar(s) || isstring(s));
addParameter(p, "AtlasPng", "", @(s) ischar(s) || isstring(s));
addParameter(p, "Transform", "affine");
addParameter(p, "Orient", "Pleft_Aright");
addParameter(p, "OutDir", "");
parse(p, varargin{:});
opt = p.Results;

root = fileparts(mfilename("fullpath"));
addpath(fileparts(root));

fovImg = ccf_dorsal_load_fov(fov, "FlipLR", opt.FlipLR);
fovShow = ccf_dorsal_enhance(fovImg);

atlasPng = char(string(opt.AtlasPng));
if strlength(strtrim(atlasPng)) > 0 && exist(atlasPng, "file") == 2
    Ref = ccf_dorsal_reference_from_png(atlasPng, opt.Orient);
else
    Ref = ccf_dorsal_make_reference("Orient", opt.Orient);
end

photoPath = char(string(opt.Photo1p));
use1p = strlength(strtrim(photoPath)) > 0;
if use1p
    photo = ccf_dorsal_load_fov(photoPath, "FlipLR", false);
    photoShow = ccf_dorsal_enhance(photo);
    fprintf("Step 1/2: click matching landmarks on MESOSCOPE FOV then 1P PHOTO.\n");
    S1 = run_cpselect(fovShow, photoShow, char(string(opt.Transform)), "fov_to_1p");
    fprintf("Step 2/2: click matching landmarks on 1P PHOTO then ALLEN DORSAL MAP.\n");
    S2 = run_cpselect(photoShow, Ref.image, char(string(opt.Transform)), "onep_to_ccf");
    S = ccf_dorsal_compose(S1, S2);
    S.photo = photoPath;
    movingShow = fovShow;
    fixedImg = Ref.image;
else
    fprintf("Click matching landmarks on FOV (left) then ALLEN DORSAL MAP (right).\n");
    fprintf("Need >= 3 points for affine, 2 for similarity, 4 for projective.\n");
    S = run_cpselect(fovShow, Ref.image, char(string(opt.Transform)), "fov_to_ccf");
    movingShow = fovShow;
    fixedImg = Ref.image;
end

S.fov_source = fov;
S.flipLR = logical(opt.FlipLR);
S.orient = char(string(opt.Orient));
S.atlas = Ref;
S.fov_size = size(fovImg);

outDir = char(string(opt.OutDir));
if strlength(strtrim(outDir)) == 0
    outDir = fullfile(fileparts(root), "ccf_dorsal_output");
end
if exist(outDir, "dir") ~= 7
    mkdir(outDir);
end
S.outDir = outDir;

Rfixed = imref2d(size(fixedImg, 1:2));
warped = imwarp(movingShow, S.tform, "OutputView", Rfixed, "SmoothEdges", true);
write_qc(movingShow, fixedImg, warped, S, outDir);
save(fullfile(outDir, "registration.mat"), "S", "-v7.3");
fprintf("Saved %s\n", fullfile(outDir, "registration.mat"));
end


function S = run_cpselect(moving, fixed, ttype, name)
if ndims(moving) < 3
    moving = repmat(moving, [1 1 3]);
end
if ndims(fixed) < 3
    fixed = repmat(fixed, [1 1 3]);
end
ccf_dorsal_show_landmark_guide();
[mp, fp] = cpselect(moving, fixed, "Wait", true);
if size(mp, 1) < 2
    error("ccf_dorsal_register:NeedPoints", "Need at least 2 paired landmarks.");
end
ttype = lower(char(ttype));
if strcmp(ttype, "projective") && size(mp, 1) < 4
    warning("ccf_dorsal_register:FewPoints", "Projective needs 4 points; using affine.");
    ttype = "affine";
end
if exist("fitgeotform2d", "file") == 2
    tform = fitgeotform2d(mp, fp, ttype);
else
    tform = fitgeotrans(mp, fp, ttype);
end
S = struct();
S.name = name;
S.moving_points = mp;
S.fixed_points = fp;
S.transform_type = ttype;
S.tform = tform;
end


function write_qc(moving, fixed, warped, S, outDir)
figure("Position", [20, 20, 1600, 620], "Color", "w", "Visible", "off");
subplot(1, 3, 1);
imshow(moving);
title("FOV (moving)");
hold on;
plot(S.moving_points(:, 1), S.moving_points(:, 2), "y+", "MarkerSize", 10, "LineWidth", 1.5);
subplot(1, 3, 2);
imshow(fixed);
title("Allen dorsal map (fixed)");
hold on;
plot(S.fixed_points(:, 1), S.fixed_points(:, 2), "c+", "MarkerSize", 10, "LineWidth", 1.5);
subplot(1, 3, 3);
imshow(imfuse(fixed, warped, "falsecolor"));
title("Overlay  (magenta = atlas, green = warped FOV)");
exportgraphics(gcf, fullfile(outDir, "qc_overlay.png"), "Resolution", 160);
close(gcf);
end
