function tblAll = ccf_label_melittin(varargin)
%CCF_LABEL_MELITTIN  Allen area for every melittin / melittin+morphine injection cell.
%
% Mirrors ccf_label_facepoke_n18 for the melittin arm. Cells come from the INJECTION
% session (baseline and injection have independent EXTRACT IDs, so there is no tracking
% across the injection - the up/down call is injection vs the baseline session mean).
%
%   tbl = ccf_label_melittin
%   tbl = ccf_label_melittin("Only", "CS1014")
%
% Writes melittin_ccf/cell_ccf_all.csv and per-recording CSVs.

p = inputParser;
p.addParameter("Only", "", @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
only = string(p.Results.Only);

here = fileparts(mfilename("fullpath"));
root = fileparts(here);
addpath(here, root);
outRoot = fullfile(here, "melittin_ccf");
if exist(outRoot, "dir") ~= 7
    mkdir(outRoot);
end

jobs = local_jobs();
if only ~= ""
    jobs = jobs(contains(jobs.mouse, only) | contains(jobs.label, only), :);
end

tblAll = table();
for i = 1:height(jobs)
    lab = jobs.label(i);
    tFile = fullfile(here, "transforms", lab + ".mat");
    if exist(tFile, "file") ~= 2
        fprintf("[%d/%d] no transform, skip  %s\n", i, height(jobs), lab);
        continue
    end
    resFile = fullfile(root, "melittin_lab_meeting", "per_recording", jobs.resLabel(i), "results.mat");
    if exist(resFile, "file") ~= 2
        fprintf("[%d/%d] no results.mat, skip  %s\n", i, height(jobs), lab);
        continue
    end
    fprintf("[%d/%d] labelling  %s\n", i, height(jobs), lab);

    T = load(tFile);
    S = load(resFile, "R");
    R = S.R;

    Ri = extract_resolve_extract_files(char(jobs.injPath(i)));
    if isempty(Ri.frr) || exist(Ri.frr, "file") ~= 2
        warning("ccf_label_melittin:NoExt", "EXTRACT missing for %s", lab);
        continue
    end

    % The melittin frr files store S as CSC under /S/{ir,jc,data}, which is not the layout
    % extract_load_ext_spatial expects. Read it the same way extract_melittin_lab_meeting
    % does, so the centroids are in exactly the frame the analysis used.
    nNeurons = numel(R.amp);
    [Sn, H] = local_frr_spatial(Ri.frr, nNeurons, 2304);
    nPix = size(Sn, 1);
    rowS = rem((0:nPix-1)', H) + 1;
    colS = floor((0:nPix-1)' / H) + 1;
    wsum = full(sum(Sn, 1))';
    cx = full((rowS' * Sn))' ./ max(wsum, eps);
    cy = full((colS' * Sn))' ./ max(wsum, eps);
    if T.flipLR
        cx = H + 1 - cx;
    end

    tbl = ccf_apply_to_cells(T, [cx, cy], "flipApplied", true);
    tbl.label = repmat(lab, height(tbl), 1);
    tbl.mouse = repmat(jobs.mouse(i), height(tbl), 1);
    tbl.cohort = repmat(jobs.cohort(i), height(tbl), 1);
    tbl.amp = R.amp(:);
    tbl.pval = R.pval(:);
    % amp is injection mean minus baseline session mean, so normalise by the baseline
    % level to make recordings with different expression comparable
    tbl.amp_pct = 100 * R.amp(:) / R.mean_dff_bl;
    tbl.is_exc = logical(R.is_up(:));
    tbl.is_inh = logical(R.is_down(:));

    recDir = fullfile(outRoot, lab);
    if exist(recDir, "dir") ~= 7
        mkdir(recDir);
    end
    writetable(tbl, fullfile(recDir, "cell_ccf.csv"));
    tblAll = [tblAll; tbl]; %#ok<AGROW>
end

if height(tblAll) == 0
    fprintf("Nothing labelled.\n");
    return
end
writetable(tblAll, fullfile(outRoot, "cell_ccf_all.csv"));
fprintf("\nwrote %s (%d cells)\n", fullfile(outRoot, "cell_ccf_all.csv"), height(tblAll));

for c = unique(tblAll.cohort)'
    sub = tblAll(tblAll.cohort == c, :);
    fprintf("  %-18s %6d cells, %d mice, up %.1f%%, down %.1f%%\n", c, height(sub), ...
        numel(unique(sub.mouse)), 100 * mean(sub.is_exc), 100 * mean(sub.is_inh));
end
end


function jobs = local_jobs()
srv = "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon";
spec = { ...
    "melittin_only",     fullfile(srv, "meltinonly"),        "CS1014-1b", "CS1014"; ...
    "melittin_only",     fullfile(srv, "meltinonly"),        "CS1014-2w", "CS1014"; ...
    "melittin_only",     srv,                                "CS0204-1b", "CS0204"; ...
    "melittin_only",     srv,                                "CS0204-2w", "CS0204"; ...
    "melittin_morphine", fullfile(srv, "morphine_20260408"), "1014m01b",  "morphine"; ...
    "melittin_morphine", fullfile(srv, "morphine_20260408"), "1014m02w",  "morphine"};

cohort = strings(0, 1); mouse = strings(0, 1); label = strings(0, 1);
resLabel = strings(0, 1); injPath = strings(0, 1);
for i = 1:size(spec, 1)
    for mir = ["mirror1", "mirror2"]
        cohort(end + 1, 1) = spec{i, 1}; %#ok<AGROW>
        mouse(end + 1, 1) = spec{i, 3}; %#ok<AGROW>
        label(end + 1, 1) = "mel_" + spec{i, 3} + "_" + mir; %#ok<AGROW>
        resLabel(end + 1, 1) = spec{i, 3} + "_" + spec{i, 4} + "_" + mir; %#ok<AGROW>
        injPath(end + 1, 1) = fullfile(spec{i, 2}, spec{i, 3}, "Injection", mir + "path"); %#ok<AGROW>
    end
end
jobs = table(label, mouse, cohort, resLabel, injPath);
end


function [Sn, H] = local_frr_spatial(frrFile, nNeurons, framePx)
ir = double(h5read(frrFile, "/S/ir"));
jc = double(h5read(frrFile, "/S/jc"));
data = double(h5read(frrFile, "/S/data"));
nCol = numel(jc) - 1;
nPix = framePx * framePx;
colIdx = zeros(jc(end), 1);
for c = 1:nCol
    colIdx(jc(c)+1:jc(c+1)) = c;
end
Ss = sparse(ir + 1, colIdx, data, nPix, nCol);
if nCol > nNeurons
    Ss = Ss(:, 1:nNeurons);
elseif nCol < nNeurons
    Ss(nPix, nNeurons) = 0;
end
cs = full(sum(Ss, 1));
cs(cs == 0) = 1;
Sn = Ss ./ cs;
H = framePx;
end
