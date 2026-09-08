function [sk_bl, ul_bl, sk_inj, ul_inj, uf_bl, uf_inj] = extract_resolve_time_window_params(params)
%EXTRACT_RESOLVE_TIME_WINDOW_PARAMS  Per-arm time crop for Baseline vs Injection traces.
%
%   skip_first_sec_*, use_last_sec_* — same as extract_apply_time_window_T (skip + last N s).
%   use_first_sec_*  — if > 0, keep FIRST N seconds only (overrides skip/last for that arm).
%
%   If any per-arm field is present, split mode. Else legacy symmetric skip_first_sec /
%   use_last_sec on both arms.

if isempty(params)
    params = struct();
end

has_split = isfield(params, 'skip_first_sec_bl') || isfield(params, 'use_last_sec_bl') ...
    || isfield(params, 'use_first_sec_bl') ...
    || isfield(params, 'skip_first_sec_inj') || isfield(params, 'use_last_sec_inj') ...
    || isfield(params, 'use_first_sec_inj');

if has_split
    sk_bl  = pget(params, 'skip_first_sec_bl', 0);
    ul_bl  = pget(params, 'use_last_sec_bl', inf);
    uf_bl  = pget(params, 'use_first_sec_bl', 0);
    sk_inj = pget(params, 'skip_first_sec_inj', 0);
    ul_inj = pget(params, 'use_last_sec_inj', inf);
    uf_inj = pget(params, 'use_first_sec_inj', 0);
else
    sk = pget(params, 'skip_first_sec', 0);
    ul = pget(params, 'use_last_sec', inf);
    sk_bl = sk;
    ul_bl = ul;
    sk_inj = sk;
    ul_inj = ul;
    uf_bl = 0;
    uf_inj = 0;
end

end

function v = pget(s, name, default)
if isstruct(s) && isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
