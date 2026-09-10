function Ref = ccf_dorsal_make_reference(varargin)
%CCF_DORSAL_MAKE_REFERENCE  Built-in top-down dorsal-cortex sketch in bregma mm.
%
% This is a landmark guide (outline + bregma/lambda/OB), not a replacement
% for AllenAtlasTopDown. For publication area names, pass a real Allen
% top-down PNG into ccf_dorsal_register(..., 'AtlasPng', path).
%
%   Ref = ccf_dorsal_make_reference('Orient', 'Pleft_Aright')
%   Ref = ccf_dorsal_make_reference('Orient', 'A_up')

p = inputParser;
addParameter(p, "Orient", "Pleft_Aright");
addParameter(p, "MmPerPx", 0.02);
parse(p, varargin{:});
orient = char(string(p.Results.Orient));
mmPerPx = double(p.Results.MmPerPx);

apMin = -4.6;
apMax = 3.8;
mlMax = 5.6;
ap = apMin:mmPerPx:apMax;
ml = -mlMax:mmPerPx:mlMax;
[ML, AP] = meshgrid(ml, ap);
% A_up image: row 1 = anterior
APup = flipud(AP);
MLup = ML;

cortex = abs(MLup) <= (4.6 + 0.35 * max(0, 2.2 - APup)) & ...
    APup <= 3.55 & APup >= -4.35 & ...
    ~(APup > 2.7 & abs(MLup) > 2.4);

img = 0.12 * ones(size(cortex));
img(cortex) = 0.92;

% midline
img(:, abs(ml) < 0.04) = min(img(:, abs(ml) < 0.04), 0.35);

% area outlines (left hemisphere, mirrored)
areas = dorsal_areas();
RGB = cat(3, img, img, img);
for i = 1:numel(areas)
    mask = poly_mask(areas(i).ml, areas(i).ap, ml, ap);
    RGB = overlay_outline(RGB, mask, areas(i).color);
end

% landmark dots in A_up pixel coords
lm = suggested_landmarks();
for i = 1:numel(lm)
    [rr, cc] = mm_to_rc(lm(i).ml, lm(i).ap, ml, ap);
    RGB = draw_dot(RGB, rr, cc, lm(i).color);
end

if strcmpi(orient, "Pleft_Aright")
    % 90 deg CW: anterior (top) -> right
    RGB = imrotate(RGB, -90);
    % rotate landmark pixels
    H = size(img, 1);
    for i = 1:numel(lm)
        [rr, cc] = mm_to_rc(lm(i).ml, lm(i).ap, ml, ap);
        rr2 = cc;
        cc2 = H + 1 - rr;
        lm(i).row = rr2;
        lm(i).col = cc2;
    end
else
    for i = 1:numel(lm)
        [rr, cc] = mm_to_rc(lm(i).ml, lm(i).ap, ml, ap);
        lm(i).row = rr;
        lm(i).col = cc;
    end
end

% labels
RGB = insert_text_safe(RGB, lm);

Ref = struct();
Ref.image = min(max(RGB, 0), 1);
Ref.mmPerPx = mmPerPx;
Ref.ml_axis = ml;
Ref.ap_axis = ap;
Ref.orient = orient;
Ref.landmarks = lm;
Ref.note = "Built-in sketch for landmark picking. Swap in an Allen top-down PNG for area IDs.";
end


function Ref = ccf_dorsal_reference_from_png(pngPath, orient)
img = im2double(imread(pngPath));
if size(img, 3) == 1
    img = repmat(img, [1 1 3]);
end
if strcmpi(char(string(orient)), "Pleft_Aright")
    % assume user PNG is A-up; rotate to match lab maps
    % if already P-left, skip by naming file *Pleft*
    if ~contains(lower(pngPath), "pleft")
        img = imrotate(img, -90);
    end
end
Ref = struct();
Ref.image = img;
Ref.mmPerPx = NaN;
Ref.ml_axis = [];
Ref.ap_axis = [];
Ref.orient = char(string(orient));
Ref.landmarks = suggested_landmarks();
Ref.note = pngPath;
end


function [rr, cc] = mm_to_rc(mlPt, apPt, ml, ap)
cc = 1 + (mlPt - ml(1)) / (ml(2) - ml(1));
rr = 1 + (ap(end) - apPt) / (ap(2) - ap(1));
rr = max(1, min(numel(ap), round(rr)));
cc = max(1, min(numel(ml), round(cc)));
end


function mask = poly_mask(mlV, apV, ml, ap)
[ML, AP] = meshgrid(ml, ap); %#ok<ASGLU>
APup = flipud(AP);
MLup = ML;
maskL = inpolygon(MLup, APup, mlV, apV);
maskR = inpolygon(MLup, APup, -mlV, apV);
mask = maskL | maskR;
end


function RGB = overlay_outline(RGB, mask, col)
if size(mask, 1) ~= size(RGB, 1) || size(mask, 2) ~= size(RGB, 2)
    return
end
e = edge(mask, "sobel");
for c = 1:3
    sl = RGB(:, :, c);
    sl(e) = col(c);
    RGB(:, :, c) = sl;
end
end


function RGB = draw_dot(RGB, rr, cc, col)
r = round(rr); c = round(cc);
[H, W, ~] = size(RGB);
[yy, xx] = ndgrid(max(1, r-6):min(H, r+6), max(1, c-6):min(W, c+6));
keep = (yy - r).^2 + (xx - c).^2 <= 5.5^2;
for k = 1:3
    sl = RGB(:, :, k);
    idx = sub2ind([H W], yy(keep), xx(keep));
    sl(idx) = col(k);
    RGB(:, :, k) = sl;
end
end


function RGB = insert_text_safe(RGB, lm)
try
    for i = 1:numel(lm)
        RGB = insertText(RGB, [lm(i).col, lm(i).row] + [8, -12], lm(i).name, ...
            "FontSize", 14, "BoxOpacity", 0, "TextColor", "black", "AnchorPoint", "LeftBottom");
    end
catch
end
end


function areas = dorsal_areas()
% Approximate dorsal patches in mm (ML negative = left, AP+ = anterior).
areas = struct("name", {}, "ml", {}, "ap", {}, "color", {});
areas(end+1) = pack_area("MOs", [-0.2 2.8 -1.4 3.4 -1.7 2.2 -0.4 1.6], [0.75 0.35 0.20]);
areas(end+1) = pack_area("MOp", [-1.4 2.6 -2.8 2.8 -2.6 1.3 -1.3 1.3], [0.85 0.55 0.20]);
areas(end+1) = pack_area("SSp", [-1.6 1.5 -3.6 1.6 -3.8 0.0 -1.8 -0.2], [0.20 0.55 0.35]);
areas(end+1) = pack_area("SSs", [-3.5 0.6 -4.6 0.4 -4.4 -0.8 -3.3 -0.6], [0.15 0.40 0.55]);
areas(end+1) = pack_area("AUD", [-4.0 -0.6 -4.8 -0.8 -4.6 -1.8 -3.6 -1.4], [0.55 0.25 0.55]);
areas(end+1) = pack_area("VIS", [-1.0 -1.6 -3.4 -1.8 -3.2 -3.6 -0.8 -3.2], [0.20 0.30 0.75]);
areas(end+1) = pack_area("RSP", [-0.15 -0.6 -1.2 -1.0 -1.1 -3.4 -0.15 -3.6], [0.55 0.20 0.25]);
areas(end+1) = pack_area("ACA", [-0.15 1.4 -0.9 2.0 -0.8 0.4 -0.15 0.3], [0.70 0.70 0.25]);
end


function a = pack_area(name, xy, col)
xy = reshape(xy, 2, [])';
a = struct("name", name, "ml", xy(:, 1)', "ap", xy(:, 2)', "color", col);
end


function lm = suggested_landmarks()
lm = struct("name", {}, "ml", {}, "ap", {}, "color", {}, "row", {}, "col", {});
lm(1) = pack_lm("L-OB", -1.3, 3.3, [0.10 0.35 0.90]);
lm(2) = pack_lm("R-OB", 1.3, 3.3, [0.10 0.35 0.90]);
lm(3) = pack_lm("bregma", 0.0, 0.0, [0.90 0.15 0.15]);
lm(4) = pack_lm("lambda", 0.0, -3.9, [0.80 0.45 0.10]);
lm(5) = pack_lm("L-barrel", -3.0, -1.4, [0.20 0.65 0.30]);
lm(6) = pack_lm("R-barrel", 3.0, -1.4, [0.20 0.65 0.30]);
end


function s = pack_lm(name, ml, ap, col)
s = struct("name", name, "ml", ml, "ap", ap, "color", col, "row", NaN, "col", NaN);
end
