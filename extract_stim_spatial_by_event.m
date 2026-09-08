function extract_stim_spatial_by_event()
%EXTRACT_STIM_SPATIAL_BY_EVENT  Peak vs pre vs post response maps.
%
% Writes stim_triggered_figures/20250813_spatial_by_event/:
%   per_recording/<label>_spatial_by_event.png
%   cohort_spatial_by_event.png

figRoot = fullfile(fileparts(mfilename("fullpath")), "stim_triggered_figures");
peakRoot = fullfile(figRoot, "20250813_lab_meeting");
preRoot = fullfile(figRoot, "20250813_lab_meeting_pre");
postRoot = fullfile(figRoot, "20250813_lab_meeting_post");
outRoot = fullfile(figRoot, "20250813_spatial_by_event");
perOut = fullfile(outRoot, "per_recording");
if exist(perOut, "dir") ~= 7
    mkdir(perOut);
end

T = readtable(fullfile(peakRoot, "lab_meeting", "recording_index.csv"), "Delimiter", ",");
n = height(T);
cmapRb = red_blue_cmap();
mapsPre = cell(n, 1);
mapsPeak = cell(n, 1);
mapsPost = cell(n, 1);
labels = cell(n, 1);

for i = 1:n
    labels{i} = T.label{i};
    mapsPre{i} = load_map(preRoot, T.movie{i});
    mapsPeak{i} = load_map(peakRoot, T.movie{i});
    mapsPost{i} = load_map(postRoot, T.movie{i});
end

for i = 1:n
    climv = shared_clim(mapsPre{i}, mapsPeak{i}, mapsPost{i});
    figure("Position", [20, 20, 1500, 480], "Color", "w");
    sgtitle(sprintf("%s  |  Spatial poke map by event  (red = excited, blue = inhibited)", labels{i}), ...
        "FontSize", 13, "FontWeight", "bold", "Interpreter", "none");
    subplot(1, 3, 1);
    plot_one(mapsPre{i}, climv, cmapRb, "pre  (onset candidate)");
    subplot(1, 3, 2);
    plot_one(mapsPeak{i}, climv, cmapRb, "peak  (brightness peak)");
    subplot(1, 3, 3);
    plot_one(mapsPost{i}, climv, cmapRb, "post  (after peak)");
    outPng = fullfile(perOut, sprintf("%s_spatial_by_event.png", labels{i}));
    exportgraphics(gcf, outPng, "Resolution", 200);
    close(gcf);
end

ncol = 3;
figure("Position", [20, 20, 1100, 2400], "Color", "w");
sgtitle({"Poke response on FOV by event", ...
    "Columns: pre | peak | post     Red = excited, blue = inhibited     FOVs not aligned"}, ...
    "FontSize", 13, "FontWeight", "bold");
for i = 1:n
    climv = shared_clim(mapsPre{i}, mapsPeak{i}, mapsPost{i});
    trio = {mapsPre{i}, mapsPeak{i}, mapsPost{i}};
    for j = 1:ncol
        subplot(n, ncol, (i - 1) * ncol + j);
        plot_one(trio{j}, climv, cmapRb, "");
        if j == 1
            ylabel(strrep(labels{i}, "_", "\_"), "FontSize", 7, "Rotation", 0, ...
                "HorizontalAlignment", "right", "VerticalAlignment", "middle");
        end
        if i == 1
            ttl = {"pre", "peak", "post"};
            title(ttl{j}, "FontSize", 10);
        end
    end
end
exportgraphics(gcf, fullfile(outRoot, "cohort_spatial_by_event.png"), "Resolution", 160);
close(gcf);
fprintf("By-event maps -> %s\n", outRoot);
end

function M = load_map(root, movie)
S = load(fullfile(root, movie, "spatial_map.mat"), "mapSmall");
M = S.mapSmall;
end

function climv = shared_clim(A, B, C)
vals = abs([A(A ~= 0); B(B ~= 0); C(C ~= 0)]);
climv = prctile(vals, 95);
if isempty(climv) || ~isfinite(climv) || climv == 0
    climv = 1;
end
end

function plot_one(M, climv, cmapRb, ttl)
imagesc(M);
colormap(gca, cmapRb);
clim([-climv, climv]);
axis image off;
if strlength(ttl) > 0
    title(ttl);
end
end

function cmapRb = red_blue_cmap()
halfN = 128;
cmapRb = [[linspace(0, 1, halfN), ones(1, halfN)]' ...
    [linspace(0, 1, halfN), linspace(1, 0, halfN)]' ...
    [ones(1, halfN), linspace(1, 0, halfN)]'];
end
