function G = grimace_load(csvFile, varargin)
%GRIMACE_LOAD  Load PainFace or GrimACE per-frame scores.
%
%   G = grimace_load(csvFile)
%   G = grimace_load(csvFile, "Source", "painface")   % 4 FAU, 0-8
%   G = grimace_load(csvFile, "Source", "grimace")    % 5 FAU, 0-10
%
% G.t_s, G.mgs, G.n_fau, G.max_fau, G.scale_max, G.source, G.file

p = inputParser;
p.addRequired("csvFile", @(s) ischar(s) || isstring(s));
p.addParameter("Source", "auto", @(s) ischar(s) || isstring(s));
p.parse(csvFile, varargin{:});
src = lower(string(p.Results.Source));
csvFile = string(csvFile);
if exist(csvFile, "file") ~= 2
    error("grimace_load:Missing", "CSV not found:\n%s", csvFile);
end

T = readtable(csvFile, "TextType", "string");
names = lower(string(T.Properties.VariableNames));
T.Properties.VariableNames = matlab.lang.makeValidName(cellstr(names));
names = string(T.Properties.VariableNames);

tS = col(T, names, ["time_s", "time", "timestamp_s", "t"]);
frame = col(T, names, ["frame", "frame_idx", "framenumber", "video_frame"]);
mgs = col(T, names, ["mgs", "score", "grimace", "grimace_score", "total"]);
nFau = col(T, names, ["n_fau", "fau_count", "nfau", "n_action_units"]);

auNames = { ...
    "orbital", ["orbital", "orbitals", "eye"]; ...
    "nose", "nose"; ...
    "cheek", ["cheek", "cheeks"]; ...
    "ear", ["ear", "ears"]; ...
    "whisker", ["whisker", "whiskers"]};
au = struct();
nAuPresent = 0;
for i = 1:size(auNames, 1)
    v = col(T, names, auNames{i, 2});
    if ~isempty(v)
        au.(auNames{i, 1}) = v;
        nAuPresent = nAuPresent + 1;
    end
end

if src == "auto"
    if nAuPresent >= 5 || any(names == "cheek")
        src = "grimace";
    else
        src = "painface";
    end
end
if src == "painface"
    maxFau = 4;
    scaleMax = 8;
else
    maxFau = 5;
    scaleMax = 10;
end

n = height(T);
if isempty(tS)
    if isempty(frame)
        tS = (0:n-1)' ;
    else
        tS = (double(frame) - 1);
    end
end
tS = double(tS(:));
if isempty(frame)
    frame = (1:n)';
else
    frame = double(frame(:));
end

if isempty(mgs) && nAuPresent > 0
    mgs = zeros(n, 1);
    f = fieldnames(au);
    nFau = zeros(n, 1);
    for i = 1:numel(f)
        v = au.(f{i});
        ok = isfinite(v) & v >= 0;
        mgs(ok) = mgs(ok) + v(ok);
        nFau = nFau + double(ok);
    end
elseif isempty(mgs)
    error("grimace_load:NoScore", ...
        "Need an mgs/score column or per-action-unit columns in:\n%s", csvFile);
else
    mgs = double(mgs(:));
end
if isempty(nFau)
    nFau = repmat(maxFau, n, 1);
else
    nFau = double(nFau(:));
end

G = struct("file", csvFile, "source", src, "max_fau", maxFau, ...
    "scale_max", scaleMax, "t_s", tS, "frame", frame, "mgs", mgs, ...
    "n_fau", nFau, "au", au);
end

function v = col(T, names, aliases)
aliases = string(aliases);
v = [];
for a = aliases
    hit = names == a;
    if any(hit)
        raw = T{:, find(hit, 1)};
        if isstring(raw) || iscellstr(raw)
            raw = str2double(string(raw));
        end
        v = double(raw);
        v(~isfinite(v) | v < 0) = NaN;
        return
    end
end
end
