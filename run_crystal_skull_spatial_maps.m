%RUN_CRYSTAL_SKULL_SPATIAL_MAPS  Paint poke amp onto FOV for 20250813 packs.
%
% Processes peak, then pre, then post lab-meeting folders. Uses Omer's
% spatial_weights CSC inside each *_ext.mat (no re-export needed).

extRoot = "\\stg-tnr50a1.stanford.edu\crystal_skull\omer_hazon\20250813";
stimMat = fullfile(extRoot, "estimated_stim_timing_20250813.mat");
figRoot = fullfile(fileparts(mfilename("fullpath")), "stim_triggered_figures");

jobs = {
    "peak", fullfile(figRoot, "20250813_lab_meeting"), "peak"
    "pre",  fullfile(figRoot, "20250813_lab_meeting_pre"), "pre"
    "post", fullfile(figRoot, "20250813_lab_meeting_post"), "post"
    };

baseParams = struct( ...
    "peri_win", [-2, 3], ...
    "base_win", [-2, 0], ...
    "resp_win", [0.4, 1.2], ...
    "alpha", 0.05, ...
    "fps_movie", 33, ...
    "frame_px", 2304);

for i = 1:size(jobs, 1)
    params = baseParams;
    params.event_field = jobs{i, 3};
    fprintf("\n########## %s ##########\n", jobs{i, 1});
    extract_stim_write_spatial_maps(jobs{i, 2}, extRoot, stimMat, params);
    close all;
end
fprintf("\nAll spatial-map jobs finished.\n");
