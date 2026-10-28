function grimace_plot(G, M, S, outDir)
%GRIMACE_PLOT  Session grimace trace, score histogram, peri-poke STA.

if nargin < 4 || strlength(string(outDir)) == 0
    outDir = pwd;
end
outDir = string(outDir);
if exist(outDir, "dir") ~= 7
    mkdir(outDir);
end

ok = isfinite(G.mgs) & G.n_fau >= G.max_fau;
fig = figure("Color", "w", "Position", [40, 40, 1200, 780], "Visible", "off");

ax1 = subplot(2, 2, [1 2]);
plot(ax1, G.t_s(ok), G.mgs(ok), "-", "Color", [0.25, 0.25, 0.25], "LineWidth", 0.6);
hold(ax1, "on");
yline(ax1, 5, "--", "Color", [0.75, 0.15, 0.15], "LineWidth", 1);
if nargin >= 3 && ~isempty(S) && isfield(S, "n_pokes")
    % poke markers drawn only if they fall in view
end
xlabel(ax1, "time (s)");
ylabel(ax1, "MGS");
title(ax1, sprintf("%s   mean MGS %.2f   high-grimace %.1f%%   MaxFAU %d/%d", ...
    G.source, M.mean_mgs, M.pct_high, M.n_max_fau, M.n_frames), ...
    "Interpreter", "none");
ylim(ax1, [0, G.scale_max]);
ax1.Box = "off";

ax2 = subplot(2, 2, 3);
if M.n_max_fau > 0
    histogram(ax2, G.mgs(ok), 0:G.scale_max, "FaceColor", [0.35, 0.45, 0.70], ...
        "EdgeColor", "none", "Normalization", "probability");
end
xlabel(ax2, "MGS");
ylabel(ax2, "fraction of MaxFAU frames");
title(ax2, sprintf("low %.0f%%  med %.0f%%  high %.0f%%", ...
    M.pct_low, M.pct_med, M.pct_high));
ax2.Box = "off";

ax3 = subplot(2, 2, 4);
if nargin >= 3 && ~isempty(S) && isfield(S, "sta") && any(isfinite(S.sta))
    plot(ax3, S.t_peri, S.sta, "-", "Color", [0.75, 0.15, 0.15], "LineWidth", 1.8);
    xline(ax3, 0, "k-");
    title(ax3, sprintf("peri-poke grimace  n=%d  Delta=%.2f  p=%.3g", ...
        S.n_used, S.mean_amp, S.p_amp));
else
    text(ax3, 0.5, 0.5, "no poke times", "HorizontalAlignment", "center");
    title(ax3, "peri-poke grimace (none for melittin BL/INJ)");
end
xlabel(ax3, "time from poke (s)");
ylabel(ax3, "MGS");
ax3.Box = "off";

sgtitle(fig, "Crystal skull  |  automated grimace (PainFace / GrimACE scores)", ...
    "FontWeight", "bold");
outFile = fullfile(outDir, "fig_grimace.png");
exportgraphics(fig, outFile, "Resolution", 150);
close(fig);
writetable(table(G.t_s, G.frame, G.mgs, G.n_fau, ...
    VariableNames=["t_s", "frame", "mgs", "n_fau"]), ...
    fullfile(outDir, "grimace_frames.csv"));
fid = fopen(fullfile(outDir, "grimace_metrics.txt"), "w");
fprintf(fid, "source        %s\n", G.source);
fprintf(fid, "file          %s\n", G.file);
fprintf(fid, "n_frames      %d\n", M.n_frames);
fprintf(fid, "n_max_fau     %d\n", M.n_max_fau);
fprintf(fid, "mean_mgs      %.4f\n", M.mean_mgs);
fprintf(fid, "pct_low       %.2f\n", M.pct_low);
fprintf(fid, "pct_med       %.2f\n", M.pct_med);
fprintf(fid, "pct_high      %.2f\n", M.pct_high);
if nargin >= 3 && ~isempty(S) && isfield(S, "n_pokes")
    fprintf(fid, "n_pokes       %d\n", S.n_pokes);
    fprintf(fid, "peri_delta    %.4f\n", S.mean_amp);
    fprintf(fid, "peri_p        %.4g\n", S.p_amp);
end
fclose(fid);
fprintf("wrote %s\n", outFile);
end
