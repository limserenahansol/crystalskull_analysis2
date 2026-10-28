function M = grimace_metrics(G, varargin)
%GRIMACE_METRICS  Mean MGS and time in low / medium / high grimace (McCoy 2024).
%
% Uses MaxFAU frames only. High grimace = MGS 5-scale_max (PainFace bins).
% Need >= MinMaxFau frames or mean_mgs is NaN.

p = inputParser;
p.addParameter("MinMaxFau", 75, @isnumeric);
p.parse(varargin{:});
minN = p.Results.MinMaxFau;

ok = isfinite(G.mgs) & G.n_fau >= G.max_fau;
mgs = G.mgs(ok);
nAll = numel(G.mgs);
nOk = numel(mgs);

low = mgs <= 2;
med = mgs >= 3 & mgs <= 4;
high = mgs >= 5;
if G.scale_max > 8
    med = mgs >= 3 & mgs <= 5;
    high = mgs >= 6;
end

M = struct();
M.source = G.source;
M.n_frames = nAll;
M.n_max_fau = nOk;
M.frac_max_fau = nOk / max(nAll, 1);
M.min_max_fau = minN;
M.enough_frames = nOk >= minN;
if M.enough_frames && nOk > 0
    M.mean_mgs = mean(mgs);
    M.pct_low = 100 * mean(low);
    M.pct_med = 100 * mean(med);
    M.pct_high = 100 * mean(high);
else
    M.mean_mgs = NaN;
    M.pct_low = NaN;
    M.pct_med = NaN;
    M.pct_high = NaN;
end
end
