function [img, meta] = ccf_facepoke_fov_image(srcFolderOrExt, recLabel)
%CCF_FACEPOKE_FOV_IMAGE  EXTRACT summary in the same frame as the n18 FOV maps.
%
%   [img, meta] = ccf_facepoke_fov_image(srcFolder)
%   [img, meta] = ccf_facepoke_fov_image(extFile, recLabel)
%
% img is the transposed summary_image (not yet fliplr). Pass meta.flipLR into
% ccf_register_fov so mirror2 matches example_fov.png (P left, A right).

root = fileparts(fileparts(mfilename("fullpath")));
addpath(root);

src = char(string(srcFolderOrExt));
meta = struct("ext_file", "", "flipLR", false, "label", "");

if nargin >= 2 && strlength(string(recLabel)) > 0
    meta.label = char(string(recLabel));
    meta.flipLR = contains(meta.label, "mirror2");
end

if isfolder(src)
    stimFile = fullfile(src, "stim_results.mat");
    if exist(stimFile, "file") ~= 2
        error("ccf_facepoke_fov_image:NoStim", "No stim_results.mat in:\n%s", src);
    end
    S = load(stimFile, "R");
    extFile = char(S.R.ext_file);
    if nargin < 2 || strlength(string(recLabel)) == 0
        [~, stem] = fileparts(src);
        meta.label = stem;
        meta.flipLR = contains(stem, "mirror2");
    end
else
    extFile = src;
end

if exist(extFile, "file") ~= 2
    error("ccf_facepoke_fov_image:NoExt", "EXTRACT file not found:\n%s", extFile);
end

img = double(h5read(extFile, "/info/summary_image")');
meta.ext_file = extFile;
meta.size = size(img);
end
