function mapSmall = extract_stim_plot_spatial_fig(P)
%EXTRACT_STIM_PLOT_SPATIAL_FIG  Paint poke amplitude onto the FOV.
%
%   mapSmall = extract_stim_plot_spatial_fig(P)
%
% P fields:
%   Sn, amp, framePx   required (Sn is nPix-by-nNeurons, columns sum to 1)
%   isExc, isInh       optional logical masks (significant cells)
%   sumimg             optional anatomy underlay
%   outFile, runId     output PNG path and title label
%
% Three panels: summary image | all-neuron signed map | significant-only map.
% Red = excited (amp > 0), blue = inhibited (amp < 0).

Sn = P.Sn;
amp = P.amp(:);
framePx = P.framePx;
nPix = framePx * framePx;
if size(Sn, 1) ~= nPix
    error("extract_stim_plot_spatial_fig:SizeMismatch", ...
        "Sn has %d rows but framePx^2 = %d.", size(Sn, 1), nPix);
end
if size(Sn, 2) ~= numel(amp)
    error("extract_stim_plot_spatial_fig:NeuronMismatch", ...
        "Sn has %d columns but amp has %d entries.", size(Sn, 2), numel(amp));
end

isExc = false(size(amp));
isInh = false(size(amp));
if isfield(P, "isExc") && ~isempty(P.isExc)
    isExc = logical(P.isExc(:));
end
if isfield(P, "isInh") && ~isempty(P.isInh)
    isInh = logical(P.isInh(:));
end
sumimg = [];
if isfield(P, "sumimg")
    sumimg = P.sumimg;
end
runId = "";
if isfield(P, "runId")
    runId = P.runId;
end
outFile = "";
if isfield(P, "outFile")
    outFile = P.outFile;
end

% CSC is stored in the un-transposed HDF5 summary orientation; this
% pipeline displays summary_image with a transpose, so match that here.
respAll = reshape(full(Sn * amp), framePx, framePx)';
ampSig = amp;
ampSig(~(isExc | isInh)) = 0;
respSig = reshape(full(Sn * ampSig), framePx, framePx)';

showPx = 768;
respAllShow = imresize(respAll, [showPx, showPx]);
respSigShow = imresize(respSig, [showPx, showPx]);
sigAbs = abs(respSigShow(respSigShow ~= 0));
if isempty(sigAbs)
    sigAbs = abs(respAllShow(respAllShow ~= 0));
end
climv = prctile(sigAbs, 95);
if isempty(climv) || ~isfinite(climv) || climv == 0
    climv = 1;
end
cmapRb = red_blue_cmap();

figure("Position", [20, 20, 1680, 560], "Color", "w");
sgtitle(sprintf("%s  |  Spatial poke-response map (red = excited, blue = inhibited)", runId), ...
    "FontSize", 13, "FontWeight", "bold", "Interpreter", "none");

subplot(1, 3, 1);
if ~isempty(sumimg)
    imagesc(imresize(double(sumimg), [showPx, showPx]));
    colormap(gca, "gray");
    axis image off;
    title("Summary image (anatomy)");
else
    axis off;
    title("No summary image");
end

subplot(1, 3, 2);
imagesc(respAllShow);
colormap(gca, cmapRb);
clim([-climv, climv]);
axis image off;
cb = colorbar;
cb.Label.String = "amp (post - pre, z)";
title(sprintf("All neurons (n=%d)", numel(amp)));

subplot(1, 3, 3);
imagesc(respSigShow);
colormap(gca, cmapRb);
clim([-climv, climv]);
axis image off;
colorbar;
title(sprintf("Significant only  exc=%d  inh=%d", sum(isExc), sum(isInh)));

if strlength(outFile) > 0
    exportgraphics(gcf, outFile, "Resolution", 200);
end

mapSmall = imresize(respAll, [256, 256]);
end

function cmapRb = red_blue_cmap()
halfN = 128;
cmapRb = [[linspace(0, 1, halfN), ones(1, halfN)]' ...
    [linspace(0, 1, halfN), linspace(1, 0, halfN)]' ...
    [ones(1, halfN), linspace(1, 0, halfN)]'];
end
