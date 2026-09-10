function jobs = ccf_n18_jobs(varargin)
%CCF_N18_JOBS  The 18 face-poke recordings (13 Aug + 18 Aug).
%
%   jobs = ccf_n18_jobs
%   jobs = ccf_n18_jobs("Only", "m1f")
%   jobs = ccf_n18_jobs("Only", "m1f_s00_pwr8.38_mirror1")

p = inputParser;
p.addParameter("Only", "");
p.parse(varargin{:});

root = fileparts(fileparts(mfilename("fullpath")));
files = {
    fullfile(root, "stim_triggered_figures", "20250813_lab_meeting", "lab_meeting", "recording_index.csv")
    fullfile(root, "stim_triggered_figures", "20250818_lab_meeting", "lab_meeting", "recording_index.csv")
    };

label = strings(0, 1);
movie = strings(0, 1);
mouse = strings(0, 1);
mirror = strings(0, 1);
srcFolder = strings(0, 1);
for i = 1:numel(files)
    if exist(files{i}, "file") ~= 2
        error("ccf_n18_jobs:NoIndex", "Missing recording index:\n%s", files{i});
    end
    T = readtable(files{i}, "Delimiter", ",");
    label = [label; string(T.label)]; %#ok<AGROW>
    movie = [movie; string(T.movie)]; %#ok<AGROW>
    mouse = [mouse; string(T.mouse)]; %#ok<AGROW>
    mirror = [mirror; string(T.mirror)]; %#ok<AGROW>
    srcFolder = [srcFolder; string(T.src_folder)]; %#ok<AGROW>
end

jobs = table(label, movie, mouse, mirror, srcFolder, ...
    VariableNames=["label", "movie", "mouse", "mirror", "src_folder"]);

only = strtrim(char(string(p.Results.Only)));
if strlength(only) == 0
    return
end
lab = jobs.label;
mou = jobs.mouse;
if any(lab == only)
    jobs = jobs(lab == only, :);
elseif any(mou == only)
    jobs = jobs(mou == only, :);
else
    jobs = jobs(contains(lab, only), :);
end
end
