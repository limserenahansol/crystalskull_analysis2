function tblAll = ccf_label_facepoke_n18(varargin)
%CCF_LABEL_FACEPOKE_N18  Assign Allen area names to every face-poke cell.
%
% Needs a ccf_transform.mat from ccf_register_facepoke_n18 for each recording.
% Writes per-recording cell_ccf.csv and a pooled summary of excited vs
% inhibited cells by area.
%
%   tbl = ccf_label_facepoke_n18
%   tbl = ccf_label_facepoke_n18("Only", "m1f")

p = inputParser;
p.addParameter("Only", "");
p.parse(varargin{:});

thisDir = fileparts(mfilename("fullpath"));
addpath(thisDir);
addpath(fileparts(thisDir));
outRoot = fullfile(thisDir, "n18");
jobs = ccf_n18_jobs("Only", p.Results.Only);

tblAll = table();
for i = 1:height(jobs)
    recLabel = char(jobs.label(i));
    tFile = fullfile(outRoot, recLabel, "ccf_transform.mat");
    if exist(tFile, "file") ~= 2
        fprintf("[%d/%d] no transform, skip  %s\n", i, height(jobs), recLabel);
        continue
    end
    fprintf("[%d/%d] labelling  %s\n", i, height(jobs), recLabel);
    T = load(tFile);
    stimFile = fullfile(char(jobs.src_folder(i)), "stim_results.mat");
    S = load(stimFile, "R");
    R = S.R;
    nNeurons = numel(R.amp);
    extFile = char(R.ext_file);
    if exist(extFile, "file") ~= 2
        warning("ccf_label_facepoke_n18:NoExt", "EXTRACT missing for %s", recLabel);
        continue
    end

    [Sn, frameSize] = extract_load_ext_spatial(extFile, nNeurons, 2304);
    H = frameSize(1);
    nPix = size(Sn, 1);
    rowS = rem((0:nPix-1)', H) + 1;
    colS = floor((0:nPix-1)' / H) + 1;
    wsum = full(sum(Sn, 1))';
    cx = full((rowS' * Sn))' ./ max(wsum, eps);
    cy = full((colS' * Sn))' ./ max(wsum, eps);
    if T.flipLR
        cx = H + 1 - cx;
    end
    cellXY = [cx, cy];

    tbl = ccf_apply_to_cells(T, cellXY, "flipApplied", true);
    tbl.label = repmat(string(recLabel), height(tbl), 1);
    tbl.mouse = repmat(jobs.mouse(i), height(tbl), 1);
    tbl.movie = repmat(jobs.movie(i), height(tbl), 1);
    tbl.amp = R.amp(:);
    tbl.pval = R.pval(:);
    tbl.is_exc = tbl.pval < 0.05 & tbl.amp > 0;
    tbl.is_inh = tbl.pval < 0.05 & tbl.amp < 0;
    writetable(tbl, fullfile(outRoot, recLabel, "cell_ccf.csv"));
    tblAll = [tblAll; tbl]; %#ok<AGROW>
end

if height(tblAll) == 0
    fprintf("Nothing labelled. Register FOVs first with ccf_register_facepoke_n18.\n");
    return
end

writetable(tblAll, fullfile(outRoot, "cell_ccf_all.csv"));
write_area_summary(tblAll, fullfile(outRoot, "area_summary.csv"));
fprintf("wrote %s (%d cells)\n", fullfile(outRoot, "cell_ccf_all.csv"), height(tblAll));
end


function write_area_summary(tbl, outFile)
fid = fopen(outFile, "w");
fprintf(fid, "mouse,kind,area,n,pct\n");
mice = unique(tbl.mouse);
kinds = {"all", "excited", "inhibited"};
for m = 1:numel(mice)
    tm = tbl(tbl.mouse == mice(m), :);
    for k = 1:numel(kinds)
        switch kinds{k}
            case "excited"
                sub = tm(tm.is_exc, :);
            case "inhibited"
                sub = tm(tm.is_inh, :);
            otherwise
                sub = tm;
        end
        if height(sub) == 0
            continue
        end
        [u, ~, g] = unique(sub.area_acronym);
        cnt = accumarray(g, 1);
        [cnt, ord] = sort(cnt, "descend");
        u = u(ord);
        for i = 1:numel(u)
            fprintf(fid, "%s,%s,%s,%d,%.3f\n", mice(m), kinds{k}, u(i), cnt(i), ...
                cnt(i) / height(sub));
        end
    end
end
fclose(fid);

fprintf("\nExcited cells by area (all registered mice):\n");
exc = tbl(tbl.is_exc, :);
if height(exc) == 0
    return
end
[u, ~, g] = unique(exc.area_acronym);
cnt = accumarray(g, 1);
[cnt, ord] = sort(cnt, "descend");
u = u(ord);
for i = 1:min(12, numel(u))
    fprintf("  %-10s %5d  (%.1f%%)\n", u(i), cnt(i), 100 * cnt(i) / height(exc));
end
end
