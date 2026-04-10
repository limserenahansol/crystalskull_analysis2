function S = extract_resolve_extract_files(folder)
%EXTRACT_RESOLVE_EXTRACT_FILES  Find EXTRACT .mat files in a session folder.
%
%   Supports:
%     *moco_frr.mat
%     *moco_ds_pre_ext.mat  or  *moco_ds_ext.mat  or  M_moco_ds_ext.mat
%     M_summary.mat  or  *...summary.mat  (excludes moco_ds / moco_frr in name)
%
%   S.frr, S.ds_ext, S.summary are full paths; '' if not found.

S = struct('frr', '', 'ds_ext', '', 'summary', '');
if exist(folder, 'dir') ~= 7
    return;
end

d = dir(fullfile(folder, '*.mat'));
names = {d.name};
if isempty(names)
    return;
end

% Full-rate traces + spatial S
idx = find(cellfun(@(n) ~isempty(regexp(n, 'moco_frr\.mat$', 'once')), names), 1);
if isempty(idx)
    return;
end
S.frr = fullfile(folder, names{idx});

% Downsampled / temporal weights
idx = find(endsWith(names, 'moco_ds_pre_ext.mat'), 1);
if isempty(idx)
    idx = find(endsWith(names, 'moco_ds_ext.mat') & ~endsWith(names, 'moco_ds_pre_ext.mat'), 1);
end
if isempty(idx)
    idx = find(strcmp(names, 'M_moco_ds_ext.mat'), 1);
end
if isempty(idx)
    S.frr = '';
    return;
end
S.ds_ext = fullfile(folder, names{idx});

% Framerate summary .mat
idx = find(strcmp(names, 'M_summary.mat'), 1);
if isempty(idx)
    for ii = 1:numel(names)
        n = names{ii};
        if numel(n) < 11
            continue;
        end
        if ~strcmpi(n(end-10:end), 'summary.mat')
            continue;
        end
        if ~isempty(strfind(n, 'moco_ds')) || ~isempty(strfind(n, 'moco_frr')) %#ok<STREMP>
            continue;
        end
        idx = ii;
        break;
    end
end
if isempty(idx)
    S.frr = '';
    S.ds_ext = '';
    return;
end
S.summary = fullfile(folder, names{idx});

end
