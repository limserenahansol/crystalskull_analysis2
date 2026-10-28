function out = ccf_consolidate_scale(varargin)
%CCF_CONSOLIDATE_SCALE  Recover the unknown um/pixel from the data, then re-solve every FOV.
%
% You do not need to know the objective's um/pixel. It is the same microscope for every
% recording, so the true scale is one number shared by all of them. Fitting it per FOV, as
% the 2-click registration does, spreads click error into the scale: the first pass gave
% 5.25-8.89 um/px, a 1.7x spread that is physically impossible.
%
% So: take the median of the per-FOV fits as the shared scale, then re-derive each
% transform with the scale held fixed, using the sinus confluence for position and the
% auto-detected midline for rotation. That throws away the noisiest degree of freedom and
% keeps the two reliable ones.
%
%   out = ccf_consolidate_scale                    % learn the scale, re-solve, report
%   out = ccf_consolidate_scale("DryRun", true)    % report only, write nothing
%   out = ccf_consolidate_scale("UmPerPx", 6.5)    % override with a known value
%   out = ccf_consolidate_scale("Only", "poke")    % restrict to a subset
%
% Originals are backed up to transforms/pre_consolidate/ the first time this runs.

p = inputParser;
p.addParameter("Only", "", @(s) ischar(s) || isstring(s));
p.addParameter("UmPerPx", NaN, @isnumeric);
p.addParameter("DryRun", false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
only = string(p.Results.Only);
dryRun = logical(p.Results.DryRun);

here = fileparts(mfilename("fullpath"));
addpath(here);
tDir = fullfile(here, "transforms");
bDir = fullfile(tDir, "pre_consolidate");
cacheDir = fullfile(here, "cache");

d = dir(fullfile(tDir, "*.mat"));
if only ~= ""
    d = d(contains(string({d.name}), only));
end
if isempty(d)
    error("ccf_consolidate_scale:None", "No transforms in %s. Run ccf_register_assisted first.", tDir);
end

A = load(fullfile(here, "allen_dorsal_map.mat"));
resUm = double(A.res_um);
lambdaRC = double(A.lambda_rc) + 1;
confluenceRC = [lambdaRC(1), round(double(A.bregma_rc(2)) + 1)];

% ---- gather the per-FOV fits
n = numel(d);
lab = strings(n, 1);
old = nan(n, 1);
S = cell(n, 1);
for k = 1:n
    S{k} = load(fullfile(d(k).folder, d(k).name));
    lab(k) = string(S{k}.label);
    if isfield(S{k}, "um_per_px")
        old(k) = S{k}.um_per_px;
    end
end

shared = p.Results.UmPerPx;
source = "supplied";
if ~isfinite(shared)
    shared = median(old, "omitnan");
    source = "median of the per-FOV fits";
end

fprintf("\n%d transform(s)\n", n);
fprintf("per-FOV fitted scale : median %.2f, range %.2f - %.2f um/px (spread %.2fx)\n", ...
    median(old, "omitnan"), min(old), max(old), max(old) / min(old));
fprintf("shared scale to use  : %.2f um/px  (%s)\n", shared, source);
fprintf("implied FOV width    : %.2f mm across %d px\n\n", shared * 2304 / 1000, 2304);

% ---- re-solve each with the scale pinned
newScale = nan(n, 1);
shiftUm = nan(n, 1);
for k = 1:n
    T = S{k};
    if ~isfield(T, "fov_points") || isempty(T.fov_points)
        continue
    end
    % the confluence is the last point clicked in either mode
    pC = T.fov_points(end, :);
    theta = deg2rad(T.midline_angle_deg);
    s = shared / resUm;
    R = s * [cos(theta), sin(theta); -sin(theta), cos(theta)];
    t = [confluenceRC(2); confluenceRC(1)] - R * pC(:);
    newT = affinetform2d([R, t; 0 0 1]);

    % how far did the FOV centre move, in atlas microns
    c = [T.fov_size(2), T.fov_size(1)] / 2;
    before = transformPointsForward(T.tform, c);
    after = transformPointsForward(newT, c);
    shiftUm(k) = norm(after - before) * resUm;
    newScale(k) = shared;

    if ~dryRun
        if exist(bDir, "dir") ~= 7
            mkdir(bDir);
        end
        bak = fullfile(bDir, d(k).name);
        if exist(bak, "file") ~= 2
            copyfile(fullfile(d(k).folder, d(k).name), bak);
        end
        T.tform = newT;
        T.um_per_px = shared;
        T.scale_source = source;
        T.scale_consolidated = true;
        T.um_per_px_before = old(k);
        save(fullfile(d(k).folder, d(k).name), "-struct", "T");
        local_qc(cacheDir, T, A, newT, fullfile(tDir, T.label + "_qc.png"));
        if startsWith(T.label, "poke_")
            recDir = fullfile(here, "n18", extractAfter(T.label, "poke_"));
            if exist(recDir, "dir") ~= 7
                mkdir(recDir);
            end
            save(fullfile(recDir, "ccf_transform.mat"), "-struct", "T");
        end
    end
end

out = table(lab, old, newScale, shiftUm, ...
    VariableNames=["label", "um_per_px_before", "um_per_px_after", "centre_shift_um"]);
out = sortrows(out, "centre_shift_um", "descend");

fprintf("%-36s %10s %10s %12s\n", "label", "before", "after", "shift (um)");
for k = 1:height(out)
    fprintf("%-36s %10.2f %10.2f %12.0f\n", out.label(k), out.um_per_px_before(k), ...
        out.um_per_px_after(k), out.centre_shift_um(k));
end
fprintf("\nmedian centre shift %.0f um, worst %.0f um\n", ...
    median(out.centre_shift_um, "omitnan"), max(out.centre_shift_um));

if dryRun
    fprintf("\nDRY RUN - nothing written.\n");
else
    fprintf("\nrewritten in %s (originals in pre_consolidate/)\n", tDir);
    fprintf("Look at the refreshed _qc.png files, then run ccf_label_facepoke_n18.\n");
end
end


function local_qc(cacheDir, T, A, tform, dst)
f = fullfile(cacheDir, T.label + ".mat");
if exist(f, "file") ~= 2
    return
end
Simg = load(f);
img = double(Simg.img);
if T.flipLR
    img = fliplr(img);
end
R = imref2d(size(img));
bw = imwarp(single(logical(A.boundary)), invert(tform), "nearest", "OutputView", R) > 0.5;
v = img(isfinite(img));
lo = prctile(v, 1); hi = prctile(v, 99.5);
g = min(max((img - lo) / max(hi - lo, 1), 0), 1);
c1 = g; c2 = g; c3 = g;
c1(bw) = 1; c2(bw) = 0.2; c3(bw) = 0.2;
fh = figure("Visible", "off", "Color", "w", "Position", [50 50 900 900]);
imshow(cat(3, c1, c2, c3));
title(sprintf("%s   %.2f um/px (shared)", T.label, T.um_per_px), "Interpreter", "none");
exportgraphics(fh, dst, "Resolution", 120);
close(fh);
end
