function extract_melittin_paired_fov_figure()
%EXTRACT_MELITTIN_PAIRED_FOV_FIGURE  CS1014 only: melittin vs morphine FOVs.

root = fileparts(mfilename("fullpath"));
per = fullfile(root, "melittin_lab_meeting", "per_recording");
outDir = fullfile(root, "mark_brief");
if exist(outDir, "dir") ~= 7
    mkdir(outDir);
end

pairs = {
    "CS1014-1b", "mirror1", "CS1014-1b_CS1014_mirror1", "1014m01b_morphine_mirror1"
    "CS1014-1b", "mirror2", "CS1014-1b_CS1014_mirror2", "1014m01b_morphine_mirror2"
    "CS1014-2w", "mirror1", "CS1014-2w_CS1014_mirror1", "1014m02w_morphine_mirror1"
    "CS1014-2w", "mirror2", "CS1014-2w_CS1014_mirror2", "1014m02w_morphine_mirror2"
    };

halfN = 128;
cmapRb = [[linspace(0, 1, halfN), ones(1, halfN)]' ...
    [linspace(0, 1, halfN), linspace(1, 0, halfN)]' ...
    [ones(1, halfN), linspace(1, 0, halfN)]'];

figure("Position", [20, 20, 1800, 820], "Color", "w", "Visible", "off");
sgtitle({"Paired CS1014 mice  |  injection FOV maps", ...
    "Red = up vs baseline session mean, blue = down. Same mouse x mirror, different day. Cells are not matched."}, ...
    "FontSize", 14, "FontWeight", "bold");

for i = 1:4
    mouse = char(pairs{i, 1});
    mir = char(pairs{i, 2});
    Sm = load(fullfile(per, char(pairs{i, 3}), "results.mat"), "R");
    Sp = load(fullfile(per, char(pairs{i, 4}), "results.mat"), "R");
    row = 1 + (i > 2);
    col0 = 1 + 2 * mod(i - 1, 2);
    maps = {Sm.R.resp_map_small, Sp.R.resp_map_small};
    dvals = [Sm.R.d_pop_dff, Sp.R.d_pop_dff];
    labs = {
        sprintf("%s  |  melittin only  |  %s", mouse, mir)
        sprintf("%s  |  + morphine  |  %s", mouse, mir)
        };
    for k = 1:2
        ax = subplot(2, 4, (row - 1) * 4 + col0 + k - 1);
        M = maps{k};
        nz = abs(M(M ~= 0));
        climv = prctile(nz, 95);
        if isempty(climv) || ~isfinite(climv) || climv == 0
            climv = 1;
        end
        imagesc(M);
        colormap(ax, cmapRb);
        clim([-climv, climv]);
        axis image off;
        title(ax, {labs{k}, sprintf("\\Delta pop dF/F = %+.2f", dvals(k))}, ...
            "FontSize", 9, "Interpreter", "tex");
    end
end

outFile = fullfile(outDir, "fig4_paired_CS1014_melittin_vs_morphine.png");
exportgraphics(gcf, outFile, "Resolution", 180);
close(gcf);
fprintf("Wrote %s\n", outFile);
end
