function Tlist = ccf_register_facepoke_n18(varargin)
%CCF_REGISTER_FACEPOKE_N18  Register each n18 FOV to the Allen dorsal CCF map.
%
% Interactive. One recording at a time. Uses the EXTRACT summary image (same
% pixels as the example_fov maps), not the 18-panel montage.
%
%   ccf_register_facepoke_n18
%   ccf_register_facepoke_n18("Only", "m1f")           % one mouse
%   ccf_register_facepoke_n18("Only", "m1f_s00_pwr8.38_mirror1")
%   ccf_register_facepoke_n18("Mode", "bregma")        % 2 clicks, faster
%   ccf_register_facepoke_n18("SkipExisting", false)   % redo
%
% For each recording MATLAB opens cpselect:
%   left  = crystal-skull FOV (P left, A right)
%   right = Allen dorsal map in the same orientation
% Click 4-6 matching points (midline sinus is the best feature), then
% File > Close. QC overlay is saved under ccf_registration/n18/<label>/.
%
% After some or all recordings are done:
%   ccf_n18_atlas_gallery
%   ccf_label_facepoke_n18

p = inputParser;
p.addParameter("Mode", "landmarks");
p.addParameter("Transform", "affine");
p.addParameter("Only", "");
p.addParameter("SkipExisting", true);
p.parse(varargin{:});
opt = p.Results;

thisDir = fileparts(mfilename("fullpath"));
addpath(thisDir);
addpath(fileparts(thisDir));

jobs = ccf_n18_jobs("Only", opt.Only);
if height(jobs) == 0
    error("ccf_register_facepoke_n18:NoJobs", "No recordings matched Only='%s'.", ...
        char(string(opt.Only)));
end

outRoot = fullfile(thisDir, "n18");
if exist(outRoot, "dir") ~= 7
    mkdir(outRoot);
end

mode = char(string(opt.Mode));
ttype = char(string(opt.Transform));
if strcmpi(mode, "bregma")
    ttype = "similarity";
end

Tlist = cell(height(jobs), 1);
nDid = 0;
nSkip = 0;
for i = 1:height(jobs)
    recLabel = char(jobs.label(i));
    recDir = fullfile(outRoot, recLabel);
    outFile = fullfile(recDir, "ccf_transform.mat");
    if logical(opt.SkipExisting) && exist(outFile, "file") == 2
        fprintf("[%d/%d] skip (exists)  %s\n", i, height(jobs), recLabel);
        Tlist{i} = load(outFile);
        nSkip = nSkip + 1;
        continue
    end
    if exist(recDir, "dir") ~= 7
        mkdir(recDir);
    end

    fprintf("\n==== [%d/%d] %s ====\n", i, height(jobs), recLabel);
    [img, meta] = ccf_facepoke_fov_image(char(jobs.src_folder(i)), recLabel);
    fprintf("  EXTRACT %s\n  flipLR=%d  size %s\n", meta.ext_file, meta.flipLR, mat2str(size(img)));

    T = ccf_register_fov(img, ...
        "mode", mode, ...
        "transform", ttype, ...
        "flipLR", meta.flipLR, ...
        "orient", "Pleft_Aright", ...
        "outFile", outFile, ...
        "saveQc", true, ...
        "closeQc", true);
    T.label = string(recLabel);
    T.mouse = jobs.mouse(i);
    T.movie = jobs.movie(i);
    T.ext_file = string(meta.ext_file);
    save(outFile, "-struct", "T");
    Tlist{i} = T;
    nDid = nDid + 1;
end

fprintf("\nRegistered %d, skipped %d, total %d.\n", nDid, nSkip, height(jobs));
fprintf("Transforms: %s\n", outRoot);
if nDid + nSkip > 0
    ccf_n18_atlas_gallery("Only", opt.Only);
end
end
