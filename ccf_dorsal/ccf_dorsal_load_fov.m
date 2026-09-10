function img = ccf_dorsal_load_fov(src, varargin)
%CCF_DORSAL_LOAD_FOV  Load a top-view image from a file, EXTRACT mat, or folder.

p = inputParser;
addParameter(p, "FlipLR", false);
parse(p, varargin{:});
flipLR = logical(p.Results.FlipLR);

if isnumeric(src)
    img = double(src);
else
    src = char(string(src));
    if exist(src, "dir") == 7
        addpath(fileparts(fileparts(mfilename("fullpath"))));
        R = extract_resolve_extract_files(src);
        if isempty(R.ds_ext)
            error("ccf_dorsal_load_fov:NoExt", "No EXTRACT ext file in:\n%s", src);
        end
        src = R.ds_ext;
    end
    [~, ~, ext] = fileparts(src);
    switch lower(ext)
        case {".png", ".tif", ".tiff", ".jpg", ".jpeg", ".bmp"}
            img = im2double(imread(src));
        case ".mat"
            try
                img = h5read(src, "/info/summary_image")';
            catch
                raw = load(src);
                img = first_numeric(raw);
            end
            img = double(img);
        otherwise
            error("ccf_dorsal_load_fov:BadType", "Unsupported input: %s", src);
    end
end

if ndims(img) == 3
    img = mean(img, 3);
end
img = double(img);
if flipLR
    img = fliplr(img);
end
end


function img = first_numeric(s)
f = fieldnames(s);
for i = 1:numel(f)
    v = s.(f{i});
    if isnumeric(v) && ndims(v) >= 2 && min(size(v)) > 8
        img = v;
        return
    end
end
error("ccf_dorsal_load_fov:NoImage", "No image array in MAT file.");
end
