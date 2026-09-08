function [Sn, frameSize, S] = extract_load_ext_spatial(extFile, nNeurons, framePx)
%EXTRACT_LOAD_EXT_SPATIAL  Reconstruct EXTRACT spatial_weights from ndSparse HDF5.
%
%   [Sn, frameSize, S] = extract_load_ext_spatial(extFile, nNeurons, framePx)
%
% Omer's *_moco_ds_pre_ext.mat files store spatial_weights as an ndSparse
% object. MATLAB cannot instantiate that class, but the CSC arrays
% (ir, jc, data) remain readable under /#refs#. This function finds the
% CSC whose column count matches nNeurons and returns:
%   S         : sparse double, nPix-by-nNeurons (raw footprint weights)
%   Sn        : same, columns normalized to sum 1 (empty columns stay 0)
%   frameSize : [height width] of the FOV
%
% nNeurons and framePx are optional. If omitted, nNeurons is taken from
% /temporal_weights and the FOV size from /info/summary_image (else framePx).

if nargin < 2
    nNeurons = [];
end
if nargin < 3
    framePx = [];
end

if exist(extFile, "file") ~= 2
    error("extract_load_ext_spatial:NoFile", "EXTRACT file not found:\n%s", extFile);
end

if isempty(nNeurons)
    twSize = double(h5info(extFile, "/temporal_weights").Dataspace.Size);
    % temporal_weights is frames x neurons (frames is the longer dim here)
    if twSize(1) > twSize(2)
        nNeurons = twSize(2);
    else
        nNeurons = twSize(1);
    end
end

frameSize = infer_frame_size(extFile, framePx);
nPix = prod(frameSize);

[row, col, val] = read_matching_csc(extFile, nNeurons, nPix);
S = sparse(row, col, val, nPix, nNeurons);

colSum = full(sum(S, 1));
Sn = S;
nonzero = colSum > 0;
if any(nonzero)
    Sn(:, nonzero) = S(:, nonzero) ./ colSum(nonzero);
end
end

function frameSize = infer_frame_size(extFile, framePx)
frameSize = [];
try
    sumSize = h5info(extFile, "/info/summary_image").Dataspace.Size;
    if numel(sumSize) >= 2 && all(sumSize(1:2) > 1)
        frameSize = double(sumSize(1:2));
    end
catch
    frameSize = [];
end
if isempty(frameSize)
    if isempty(framePx)
        framePx = 2304;
    end
    frameSize = [framePx, framePx];
end
end

function [row, col, val] = read_matching_csc(extFile, nNeurons, nPix)
refs = h5info(extFile, "/#refs#");
bestNnz = -1;
bestPath = '';
for i = 1:numel(refs.Groups)
    group = refs.Groups(i);
    names = {group.Datasets.Name};
    if ~all(ismember({'ir', 'jc', 'data'}, names))
        continue
    end
    jcSize = dataset_numel(group, "jc");
    if jcSize ~= nNeurons + 1
        continue
    end
    nnzHere = dataset_numel(group, "data");
    if nnzHere > bestNnz
        bestNnz = nnzHere;
        bestPath = group.Name;
    end
end
if bestNnz < 0
    error("extract_load_ext_spatial:NoCsc", ...
        ["No /#refs# CSC with %d columns was found in:\n%s\n" ...
        "spatial_weights may be missing or stored differently."], ...
        nNeurons, extFile);
end

ir = double(h5read(extFile, [bestPath '/ir']));
jc = double(h5read(extFile, [bestPath '/jc']));
val = double(h5read(extFile, [bestPath '/data']));
ir = ir(:);
jc = jc(:);
val = val(:);

if min(ir) == 0
    row = ir + 1;
else
    row = ir;
end
if any(row < 1) || any(row > nPix)
    error("extract_load_ext_spatial:RowRange", ...
        "CSC row indices [%g %g] fall outside 1:%d pixels.", min(row), max(row), nPix);
end

col = zeros(numel(ir), 1);
for c = 1:nNeurons
    a = jc(c) + 1;
    b = jc(c + 1);
    if b >= a
        col(a:b) = c;
    end
end
keep = col > 0;
row = row(keep);
col = col(keep);
val = val(keep);
end

function n = dataset_numel(group, name)
n = 0;
for i = 1:numel(group.Datasets)
    if strcmp(group.Datasets(i).Name, name)
        n = prod(double(group.Datasets(i).Dataspace.Size));
        return
    end
end
end
