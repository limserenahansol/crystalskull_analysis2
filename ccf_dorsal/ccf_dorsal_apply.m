function [warped, xyCcf] = ccf_dorsal_apply(imgOrXy, S, varargin)
%CCF_DORSAL_APPLY  Warp an image, or map neuron XY (pixels) into atlas pixels.
%
%   warped = ccf_dorsal_apply(fovImg, S)
%   [~, xyAtlas] = ccf_dorsal_apply(xyPx, S)   % xyPx is N-by-2 [x y] in FOV

p = inputParser;
addParameter(p, "OutputView", []);
parse(p, varargin{:});

isPts = ismatrix(imgOrXy) && size(imgOrXy, 2) == 2 && size(imgOrXy, 1) < 20000;
if isPts
    xyCcf = transformPointsForward(S.tform, imgOrXy);
    warped = [];
    return
end

img = imgOrXy;
if isempty(p.Results.OutputView)
    if isfield(S, "atlas") && isfield(S.atlas, "image")
        sz = size(S.atlas.image);
        R = imref2d(sz(1:2));
    else
        sz = size(img);
        R = imref2d(sz(1:2));
    end
else
    R = p.Results.OutputView;
end
warped = imwarp(img, S.tform, "OutputView", R, "SmoothEdges", true);
xyCcf = [];
end
