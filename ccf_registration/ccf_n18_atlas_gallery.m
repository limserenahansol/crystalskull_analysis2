function ccf_n18_atlas_gallery(varargin)
%CCF_N18_ATLAS_GALLERY  Warp registered n18 FOVs onto the Allen dorsal map.
%
%   ccf_n18_atlas_gallery
%   ccf_n18_atlas_gallery("Only", "m1f")
%
% Writes:
%   ccf_registration/n18/ccf_n18_gallery.png           18-panel, A-up
%   ccf_registration/n18/ccf_n18_overlay_by_mouse.png  one atlas per mouse

p = inputParser;
p.addParameter("Only", "");
p.parse(varargin{:});

thisDir = fileparts(mfilename("fullpath"));
addpath(thisDir);
addpath(fileparts(thisDir));
outRoot = fullfile(thisDir, "n18");
jobs = ccf_n18_jobs("Only", p.Results.Only);

ready = false(height(jobs), 1);
for i = 1:height(jobs)
    ready(i) = exist(fullfile(outRoot, char(jobs.label(i)), "ccf_transform.mat"), "file") == 2;
end
jobs = jobs(ready, :);
if height(jobs) == 0
    fprintf("No ccf_transform.mat files yet. Run ccf_register_facepoke_n18 first.\n");
    return
end

A = load(fullfile(thisDir, "allen_dorsal_map.mat"));
atlasRgb = local_atlas_rgb(A);
Rfixed = imref2d(size(A.label_img));

n = height(jobs);
ncol = min(6, n);
nrow = ceil(n / ncol);
fig = figure("Color", "w", "Position", [20 20 300*ncol 280*nrow], "Visible", "off");
sgtitle("Face-poke FOVs warped onto Allen CCF dorsal surface (A-up)", "FontSize", 14);

for i = 1:n
    recLabel = char(jobs.label(i));
    T = load(fullfile(outRoot, recLabel, "ccf_transform.mat"));
    [img, ~] = ccf_facepoke_fov_image(char(jobs.src_folder(i)), recLabel);
    fov = img;
    if T.flipLR
        fov = fliplr(fov);
    end
    fov = local_stretch(fov);
    warped = imwarp(fov, T.tform, "OutputView", Rfixed, "SmoothEdges", true);
    ax = subplot(nrow, ncol, i);
    imshow(imfuse(atlasRgb, warped, "falsecolor"));
    title(ax, recLabel, "FontSize", 8, "Interpreter", "none");
    imwrite(imfuse(atlasRgb, warped, "falsecolor"), ...
        fullfile(outRoot, recLabel, "fov_on_atlas.png"));
end
exportgraphics(fig, fullfile(outRoot, "ccf_n18_gallery.png"), "Resolution", 140);
close(fig);
fprintf("gallery -> %s\n", fullfile(outRoot, "ccf_n18_gallery.png"));

mice = unique(jobs.mouse);
figM = figure("Color", "w", "Position", [20 20 420*numel(mice) 520], "Visible", "off");
for m = 1:numel(mice)
    acc = zeros(size(A.label_img));
    hit = jobs.mouse == mice(m);
    idx = find(hit);
    for k = 1:numel(idx)
        recLabel = char(jobs.label(idx(k)));
        T = load(fullfile(outRoot, recLabel, "ccf_transform.mat"));
        [img, ~] = ccf_facepoke_fov_image(char(jobs.src_folder(idx(k))), recLabel);
        fov = img;
        if T.flipLR
            fov = fliplr(fov);
        end
        fov = local_stretch(fov);
        warped = imwarp(fov, T.tform, "OutputView", Rfixed, "SmoothEdges", true);
        acc = acc + warped;
    end
    if max(acc(:)) > 0
        acc = acc / max(acc(:));
    end
    ax = subplot(1, numel(mice), m);
    imshow(imfuse(atlasRgb, acc, "falsecolor"));
    title(ax, sprintf("%s  (%d FOVs)", mice(m), numel(idx)), "FontSize", 12, "Interpreter", "none");
end
sgtitle("Per-mouse FOVs on Allen CCF (magenta = atlas, green = warped FOV)", "FontSize", 13);
exportgraphics(figM, fullfile(outRoot, "ccf_n18_overlay_by_mouse.png"), "Resolution", 140);
close(figM);
end


function bd = local_atlas_rgb(A)
d = double(A.surface_depth_um);
m = isfinite(d) & d > 0;
g = zeros(size(d));
if any(m(:))
    g(m) = 1 - local_stretch(d(m));
end
bd = repmat(g, 1, 1, 3);
b = logical(A.boundary);
ch1 = bd(:, :, 1); ch2 = bd(:, :, 2); ch3 = bd(:, :, 3);
ch1(b) = 1; ch2(b) = 0.2; ch3(b) = 0.2;
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
out = (img - lo) / (hi - lo);
out(~isfinite(out)) = 0;
out = min(max(out, 0), 1);
out = out .^ 0.5;
end
