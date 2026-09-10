function test_ccf_roundtrip()
%TEST_CCF_ROUNDTRIP  Synthetic check that the registration recovers a known transform.
%
% Builds a fake FOV by pushing the atlas through a known similarity transform, hands
% ccf_register_fov a handful of matched points, and checks that (a) the recovered transform
% matches the truth and (b) cells land in the areas they were drawn from.

root = fileparts(mfilename("fullpath"));
addpath(root);
A = load(fullfile(root, "allen_dorsal_map.mat"));
resUm = double(A.res_um);

% ---- ground-truth transform: atlas pixels -> fake FOV pixels
theta = deg2rad(7);
scale = 2.3;                       % FOV px per atlas px
tx = 120; ty = -60;
Rm = scale * [cos(theta), -sin(theta); sin(theta), cos(theta)];
truth = affinetform2d([Rm, [tx; ty]; 0 0 1]);

% ---- a fake FOV image, just so the QC path runs on something real
depth = double(A.surface_depth_um);
depth(~isfinite(depth)) = 0;
fov = imwarp(depth, truth, "linear", "OutputView", imref2d([1400 1400]));

% ---- pick landmark areas, take their atlas centroids as the "clicked" atlas points
P = ccf_area_centroid(["VISp", "SSp-bfd", "MOp", "RSPd"], "mapFile", ...
    fullfile(root, "allen_dorsal_map.mat"));
ccfPts = [P.x, P.y];
fovPts = transformPointsForward(truth, ccfPts);     % where they would appear in the FOV

fprintf("=== round trip: %d landmark points ===\n", size(ccfPts, 1));

T = ccf_register_fov(fov, "mode", "points", "transform", "similarity", ...
    "fovPoints", fovPts, "ccfPoints", ccfPts, ...
    "mapFile", fullfile(root, "allen_dorsal_map.mat"), ...
    "outFile", fullfile(root, "test_transform.mat"), ...
    "closeQc", true);

% ---- (a) does the recovered transform match the inverse of the truth?
probe = [200 200; 900 400; 600 1100; 1200 900];
gotAtlas = transformPointsForward(T.tform, probe);
wantAtlas = transformPointsForward(invert(truth), probe);
errUm = sqrt(sum((gotAtlas - wantAtlas).^2, 2)) * resUm;
fprintf("\ntransform recovery error at 4 probe points (um): %s\n", ...
    strjoin(compose("%.2f", errUm), ", "));
assert(max(errUm) < 1, "transform not recovered: max error %.2f um", max(errUm));
fprintf("PASS  transform recovered to < 1 um\n");

% ---- (b) do cells land back in the area they came from?
rng(0);
truthAcr = strings(0, 1);
cellXY = zeros(0, 2);
for k = 1:height(P)
    m = double(A.label_img) == double(A.ids(find(string(A.acronyms(:)) == P.acronym(k), 1)));
    if P.hemisphere(k) == "L"
        m(:, ceil(double(A.bregma_rc(2)) + 1):end) = false;
    else
        m(:, 1:floor(double(A.bregma_rc(2)) + 1)) = false;
    end
    [rr, cc] = find(m);
    pick = randsample(numel(rr), min(40, numel(rr)));
    xyAtlas = [cc(pick), rr(pick)];
    cellXY = [cellXY; transformPointsForward(truth, xyAtlas)]; %#ok<AGROW>
    truthAcr = [truthAcr; repmat(P.acronym(k), numel(pick), 1)]; %#ok<AGROW>
end

tbl = ccf_apply_to_cells(T, cellXY, "mapFile", fullfile(root, "allen_dorsal_map.mat"));
agree = mean(tbl.area_acronym == truthAcr);
fprintf("\ncell area labels agreeing with ground truth: %.1f%% of %d cells\n", ...
    100 * agree, numel(truthAcr));
assert(agree > 0.98, "area labels disagree too often: %.1f%%", 100 * agree);
fprintf("PASS  area labelling round-trips\n");

% ---- sanity on the stereotaxic coordinates
v1 = tbl(tbl.area_acronym == "VISp", :);
fprintf("\nVISp cells: AP %.2f to %.2f mm, ML %.2f to %.2f mm\n", ...
    min(v1.AP_mm), max(v1.AP_mm), min(v1.ML_mm), max(v1.ML_mm));
assert(all(v1.AP_mm < 0), "VISp should be posterior to bregma (AP negative)");
fprintf("PASS  VISp is posterior to bregma\n");

delete(fullfile(root, "test_transform.mat"));
qcPng = fullfile(root, "test_transform_qc.png");
if exist(qcPng, "file") == 2
    delete(qcPng);
end
fprintf("\nALL CHECKS PASSED\n");
end
