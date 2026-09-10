function img = ccf_dorsal_enhance(img)
%CCF_DORSAL_ENHANCE  Contrast stretch for landmark picking (same idea as FOV maps).

img = double(img);
if ndims(img) == 3
    img = mean(img, 3);
end
vals = img(isfinite(img));
if isempty(vals)
    img = zeros(size(img));
    return
end
lo = prctile(vals, 2);
hi = prctile(vals, 97);
if ~(isfinite(lo) && isfinite(hi)) || hi <= lo
    lo = min(vals);
    hi = max(vals);
end
img = (img - lo) / max(hi - lo, eps);
img = min(max(img, 0), 1);
img = img .^ 0.5;
try
    img = adapthisteq(img, "NumTiles", [12, 12], "ClipLimit", 0.03, "Range", "original");
catch
end
img = min(max(img, 0), 1);
end
