function T = extract_apply_time_window_T(T, fs, skip_first_sec, use_last_sec)
%EXTRACT_APPLY_TIME_WINDOW_T  Crop fluorescence matrix T (frames × neurons) in time.
%
%   T = extract_apply_time_window_T(T, fs, skip_first_sec, use_last_sec)
%
%   skip_first_sec  — drop this many seconds from the start (0 = none).
%   use_last_sec    — after that, keep only the last this many seconds of the
%                     remainder. Inf or <=0 means keep the entire remainder.
%
%   Example (symmetric crop on both traces in older scripts): skip=600, use_last=600.
%
%   Example (injection only, last 10 min): skip_first_sec_inj=0, use_last_sec_inj=600
%            in extract_resolve_time_window_params / extract_analysis2_core — baseline
%            stays full; injection uses this helper with skip=0, use_last=600.

[nf, ~] = size(T);
fs = double(fs);
if numel(fs) > 1
    fs = fs(1);
end

skip_first_sec = double(skip_first_sec);
if skip_first_sec < 0
    skip_first_sec = 0;
end

skip_f = round(skip_first_sec * fs);
i0 = skip_f + 1;
if i0 > nf
    error('extract_apply_time_window_T:SkipTooLong', ...
        'skip_first_sec (%.1f s → %d frames) exceeds trace length (%d frames).', ...
        skip_first_sec, skip_f, nf);
end

if ~isfinite(use_last_sec) || use_last_sec <= 0
    T = T(i0:end, :);
    return;
end

use_last_sec = double(use_last_sec);
rem_len = nf - i0 + 1;
want_f = max(1, round(use_last_sec * fs));
take = min(want_f, rem_len);
i1 = nf - take + 1;
T = T(i1:end, :);

end
