function [sk_bl, ul_bl, sk_inj, ul_inj] = extract_resolve_time_window_params(params)
%EXTRACT_RESOLVE_TIME_WINDOW_PARAMS  Per-arm time crop for Baseline vs Injection traces.
%
%   If any of skip_first_sec_bl, use_last_sec_bl, skip_first_sec_inj, use_last_sec_inj
%   is present, those are used (missing arms default to skip=0, use_last=inf = full).
%
%   Otherwise legacy symmetric skip_first_sec / use_last_sec apply to BOTH arms.

if isempty(params)
    params = struct();
end

has_split = isfield(params, 'skip_first_sec_bl') || isfield(params, 'use_last_sec_bl') ...
    || isfield(params, 'skip_first_sec_inj') || isfield(params, 'use_last_sec_inj');

if has_split
    sk_bl  = pget(params, 'skip_first_sec_bl', 0);
    ul_bl  = pget(params, 'use_last_sec_bl', inf);
    sk_inj = pget(params, 'skip_first_sec_inj', 0);
    ul_inj = pget(params, 'use_last_sec_inj', inf);
else
    sk = pget(params, 'skip_first_sec', 0);
    ul = pget(params, 'use_last_sec', inf);
    sk_bl = sk;
    ul_bl = ul;
    sk_inj = sk;
    ul_inj = ul;
end

end

function v = pget(s, name, default)
if isstruct(s) && isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
