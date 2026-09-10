function ccf_cache_fov_images()
%CCF_CACHE_FOV_IMAGES  Pull every FOV summary image once and work out how many are distinct.
%
% Registration is the only manual step in this pipeline, so the useful question is not
% "how do I click 30 times" but "how many of these 30 are actually the same field of view".
% Recordings from one mouse on one day, or the two mirrors of one movie, usually share a FOV.
% Where they do, one manual registration propagates to the rest by FOV-to-FOV alignment,
% which is same-modality and therefore reliable.
%
% Writes cache/<label>.mat (+ .png) and cache/fov_groups.csv.

root = fileparts(fileparts(mfilename("fullpath")));
here = fileparts(mfilename("fullpath"));
addpath(root, here);
outDir = fullfile(here, "cache");
if exist(outDir, "dir") ~= 7
    mkdir(outDir);
end

recs = local_all_recordings(root);
fprintf("%d recordings to cache\n\n", numel(recs));

imgs = cell(numel(recs), 1);
ok = false(numel(recs), 1);
for k = 1:numel(recs)
    r = recs(k);
    dst = fullfile(outDir, r.label + ".mat");
    try
        if exist(dst, "file") == 2
            S = load(dst);
            img = S.img;
        else
            img = double(h5read(r.ext_file, "/info/summary_image")');
            save(dst, "img", "-v7");
            imwrite(local_stretch(img), fullfile(outDir, r.label + ".png"));
        end
        % mirror2 is flipped in every downstream script; do the same here so the
        % comparison is between images in the frame we will actually register
        if r.flipLR
            img = fliplr(img);
        end
        imgs{k} = img;
        ok(k) = true;
        fprintf("  [%2d/%2d] %-42s %s\n", k, numel(recs), r.label, mat2str(size(img)));
    catch ME
        fprintf(2, "  [%2d/%2d] %-42s FAILED: %s\n", k, numel(recs), r.label, ME.message);
    end
end

recs = recs(ok);
imgs = imgs(ok);
n = numel(imgs);
fprintf("\ncached %d images\n", n);

% ---- pairwise similarity on heavily downsampled, normalised images
small = cell(n, 1);
for k = 1:n
    s = imresize(local_stretch(imgs{k}), [128 128]);
    s = s - mean(s(:));
    small{k} = s / max(norm(s(:)), eps);
end
S = eye(n);
for i = 1:n
    for j = i + 1:n
        S(i, j) = sum(small{i}(:) .* small{j}(:));
        S(j, i) = S(i, j);
    end
end

% ---- group: same FOV if correlation is high
thr = 0.80;
grp = zeros(n, 1);
g = 0;
for i = 1:n
    if grp(i) > 0
        continue
    end
    g = g + 1;
    stack = i;
    while ~isempty(stack)
        a = stack(end);
        stack(end) = [];
        if grp(a) > 0
            continue
        end
        grp(a) = g;
        stack = [stack, find(grp == 0 & S(a, :)' > thr)']; %#ok<AGROW>
    end
end

fprintf("\n%d distinct FOV group(s) at r > %.2f\n", g, thr);
labels = arrayfun(@(r) r.label, recs);
cohorts = arrayfun(@(r) r.cohort, recs);
for k = 1:g
    idx = find(grp == k);
    [~, best] = max(sum(S(idx, idx), 2));      % most representative member
    fprintf("\n  group %d  (%d recordings)  reference: %s\n", k, numel(idx), labels(idx(best)));
    for i = idx'
        fprintf("      %-42s %s\n", labels(i), cohorts(i));
    end
end

isRef = false(n, 1);
for k = 1:g
    idx = find(grp == k);
    [~, best] = max(sum(S(idx, idx), 2));
    isRef(idx(best)) = true;
end

T = table(labels(:), cohorts(:), grp, isRef, ...
    VariableNames=["label", "cohort", "fov_group", "is_reference"]);
writetable(T, fullfile(outDir, "fov_groups.csv"));
save(fullfile(outDir, "fov_similarity.mat"), "S", "labels", "grp", "isRef", "thr");

fprintf("\n-> %d manual registration(s) needed, not %d.\n", g, n);
fprintf("   Register the reference of each group, then propagate with ccf_propagate_transform.\n");
fprintf("   wrote %s\n", fullfile(outDir, "fov_groups.csv"));

% ---- contact sheet
f = figure("Color", "w", "Position", [40 40 1700 950], "Visible", "off");
tl = tiledlayout(f, ceil(n / 6), 6, "TileSpacing", "compact", "Padding", "compact");
for k = 1:n
    ax = nexttile(tl);
    imshow(local_stretch(imgs{k}), "Parent", ax);
    title(ax, sprintf("g%d %s", grp(k), labels(k)), "FontSize", 6.5, "Interpreter", "none");
    if isRef(k)
        ax.XColor = "r"; ax.YColor = "r"; ax.LineWidth = 2; axis(ax, "on");
        set(ax, "XTick", [], "YTick", []);
    end
end
title(tl, sprintf("All FOVs, grouped (red = group reference, the one to register by hand) - %d groups", g));
exportgraphics(f, fullfile(outDir, "fov_contact_sheet.png"), "Resolution", 130);
close(f);
fprintf("   wrote %s\n", fullfile(outDir, "fov_contact_sheet.png"));
end


% ====================================================================== helpers
function recs = local_all_recordings(root)
recs = struct("label", {}, "ext_file", {}, "cohort", {}, "flipLR", {});

% ---- face poke, from the two recording indices
jobs = ccf_n18_jobs();
for k = 1:height(jobs)
    src = char(jobs.src_folder(k));
    stimFile = fullfile(src, "stim_results.mat");
    if exist(stimFile, "file") ~= 2
        fprintf(2, "  skip %s (no stim_results.mat)\n", jobs.label(k));
        continue
    end
    S = load(stimFile, "R");
    recs(end + 1) = struct("label", "poke_" + jobs.label(k), ...
        "ext_file", string(S.R.ext_file), "cohort", "face_poke", ...
        "flipLR", contains(jobs.label(k), "mirror2")); %#ok<AGROW>
end

% ---- melittin, injection movies (the cells we label come from the injection session)
mel = { ...
    "melittin_only", "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\meltinonly", "CS1014-1b", "CS1014"; ...
    "melittin_only", "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\meltinonly", "CS1014-2w", "CS1014"; ...
    "melittin_only", "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon", "CS0204-1b", "CS0204"; ...
    "melittin_only", "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon", "CS0204-2w", "CS0204"; ...
    "melittin_morphine", "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\morphine_20260408", "1014m01b", "morphine"; ...
    "melittin_morphine", "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\morphine_20260408", "1014m02w", "morphine"};
for i = 1:size(mel, 1)
    for mirror = ["mirror1path", "mirror2path"]
        injPath = fullfile(mel{i, 2}, mel{i, 3}, "Injection", mirror);
        try
            Ri = extract_resolve_extract_files(injPath);
        catch
            continue
        end
        if isempty(Ri.ds_ext) || exist(Ri.ds_ext, "file") ~= 2
            continue
        end
        mirShort = extractBefore(mirror, "path");
        recs(end + 1) = struct( ...
            "label", "mel_" + mel{i, 3} + "_" + mirShort, ...
            "ext_file", string(Ri.ds_ext), "cohort", string(mel{i, 1}), ...
            "flipLR", mirror == "mirror2path"); %#ok<AGROW>
    end
end
end


function out = local_stretch(img)
v = img(isfinite(img));
if isempty(v)
    out = zeros(size(img));
    return
end
lo = prctile(v, 1);
hi = prctile(v, 99.5);
if hi <= lo
    hi = lo + 1;
end
out = min(max((img - lo) / (hi - lo), 0), 1);
out(~isfinite(out)) = 0;
end
