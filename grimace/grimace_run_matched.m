function grimace_run_matched(scoreDir)
%GRIMACE_RUN_MATCHED  Load GrimACE CSVs and plot BL vs INJ per mouse/date.
%
%   grimace_run_matched
%   grimace_run_matched("grimace/grimace_output")

if nargin < 1 || strlength(string(scoreDir)) == 0
    scoreDir = fullfile(fileparts(mfilename("fullpath")), "grimace_output");
end
scoreDir = string(scoreDir);
matchFile = fullfile(fileparts(mfilename("fullpath")), "behavior_calcium_match.csv");
allFile = fullfile(scoreDir, "grimace_all_frames.csv");
if exist(allFile, "file") ~= 2
    error("grimace_run_matched:Missing", ...
        "No %s. Run python grimace/score_grimace_videos.py first.", allFile);
end

addpath(fileparts(mfilename("fullpath")));
G = grimace_load(allFile, "Source", "grimace");
M = grimace_metrics(G, "MinMaxFau", 20);
outDir = fullfile(scoreDir, "matlab_plots");
if ~exist(outDir, "dir")
    mkdir(outDir);
end
grimace_plot(G, M, struct(), outDir);

T = readtable(allFile, "TextType", "string");
if exist(matchFile, "file") == 2
    matchT = readtable(matchFile, "TextType", "string");
    writetable(matchT, fullfile(outDir, "behavior_calcium_match_copy.csv"));
end

if ismember("epoch", T.Properties.VariableNames) && ismember("mouse_id", T.Properties.VariableNames)
    local_bl_inj_plot(T, outDir);
end
fprintf("grimace_run_matched: mean MGS %.2f, MaxFAU %d/%d -> %s\n", ...
    M.mean_mgs, M.n_max_fau, M.n_frames, outDir);
end


function local_bl_inj_plot(T, outDir)
ok = ismember(lower(T.epoch), ["baseline", "injection"]) & isfinite(T.mgs) & T.n_fau >= 5;
S = T(ok, :);
if height(S) < 4
    return
end
keys = S.behavior_date + " " + S.mouse_id;
uk = unique(keys);
fig = figure("Color", "w", "Position", [80 80 900 420]);
hold on
xticks = [];
xtl = strings(0, 1);
for i = 1:numel(uk)
    sub = S(keys == uk(i), :);
    bl = sub.mgs(lower(sub.epoch) == "baseline");
    inj = sub.mgs(lower(sub.epoch) == "injection");
    x0 = (i - 1) * 3;
    if ~isempty(bl)
        scatter(x0 + 0.2 * randn(numel(bl), 1), bl, 18, "b", "filled", "MarkerFaceAlpha", 0.35);
        plot(x0, mean(bl), "bd", "MarkerFaceColor", "b", "MarkerSize", 8);
    end
    if ~isempty(inj)
        scatter(x0 + 1 + 0.2 * randn(numel(inj), 1), inj, 18, "r", "filled", "MarkerFaceAlpha", 0.35);
        plot(x0 + 1, mean(inj), "rd", "MarkerFaceColor", "r", "MarkerSize", 8);
    end
    xticks(end + 1) = x0 + 0.5; %#ok<AGROW>
    xtl(end + 1, 1) = uk(i); %#ok<AGROW>
end
xlim([-1, numel(uk) * 3 - 1])
ylabel("GrimACE MGS (sum of 5 AUs, 0-10)")
title("Sampled GrimACE: Baseline vs Injection (provisional; ears often not rateable)")
set(gca, "XTick", xticks, "XTickLabel", xtl, "XTickLabelRotation", 25)
legend({"BL frames", "BL mean", "INJ frames", "INJ mean"}, "Location", "best")
exportgraphics(fig, fullfile(outDir, "fig_grimace_bl_vs_inj.png"), "Resolution", 150);
close(fig)
end
