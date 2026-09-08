function ml = extract_cohort_normalize_mirror_list(mirror_spec)
%EXTRACT_COHORT_NORMALIZE_MIRROR_LIST  Cell list of focal-plane folder names.
%   mirror_spec — char, string, or cellstr of names e.g. 'mirror1path' or {'m1','m2'}.

if ischar(mirror_spec)
    ml = {mirror_spec};
elseif isstring(mirror_spec)
    ml = cellstr(mirror_spec);
elseif iscell(mirror_spec)
    ml = mirror_spec(:)';
else
    error('extract_cohort_normalize_mirror_list:BadType', ...
        'mirror_spec must be char, string, or cell array of names.');
end

end
