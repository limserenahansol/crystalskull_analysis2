function S = grimace_align_pokes(G, pokeTimesS, varargin)
%GRIMACE_ALIGN_POKES  Peri-poke grimace, same windows as the neural STA.
%
%   S = grimace_align_pokes(G, pokeTimesS)
%   S = grimace_align_pokes(G, pokeTimesS, "PeriWin", [-2 3], ...
%           "BaseWin", [-2 0], "RespWin", [0.4 1.2])
%
% pokeTimesS is seconds from video t=0 (same clock as G.t_s).

p = inputParser;
p.addParameter("PeriWin", [-2, 3], @isnumeric);
p.addParameter("BaseWin", [-2, 0], @isnumeric);
p.addParameter("RespWin", [0.4, 1.2], @isnumeric);
p.parse(varargin{:});
peri = p.Results.PeriWin;
base = p.Results.BaseWin;
resp = p.Results.RespWin;

ok = isfinite(G.mgs) & G.n_fau >= G.max_fau;
t = G.t_s(ok);
y = G.mgs(ok);
pokes = double(pokeTimesS(:));
pokes = pokes(isfinite(pokes));

dt = median(diff(t));
if ~(isfinite(dt) && dt > 0)
    dt = 1;
end
tGrid = peri(1):dt:peri(2);
nT = numel(tGrid);
nP = numel(pokes);
trials = nan(nP, nT);
amp = nan(nP, 1);
for k = 1:nP
    tk = t - pokes(k);
    in = tk >= peri(1) & tk <= peri(2);
    if nnz(in) < 3
        continue
    end
    trials(k, :) = interp1(tk(in), y(in), tGrid, "linear", NaN);
    b = tk >= base(1) & tk < base(2);
    r = tk >= resp(1) & tk < resp(2);
    if any(b) && any(r)
        amp(k) = mean(y(r), "omitnan") - mean(y(b), "omitnan");
    end
end

S = struct();
S.t_peri = tGrid(:);
S.n_pokes = nP;
S.n_used = sum(any(isfinite(trials), 2));
S.sta = mean(trials, 1, "omitnan")';
S.amp = amp;
S.mean_amp = mean(amp, "omitnan");
if sum(isfinite(amp)) >= 6
    S.p_amp = signrank(amp(isfinite(amp)));
else
    S.p_amp = NaN;
end
S.peri_win = peri;
S.base_win = base;
S.resp_win = resp;
end
