function tbl = ccf_apply_to_cells(T, cellXY, varargin)
%CCF_APPLY_TO_CELLS  Label EXTRACT cells with Allen CCF areas using a fitted transform.
%
%   tbl = ccf_apply_to_cells(T, cellXY)
%   tbl = ccf_apply_to_cells("ccf_transform_20260910_141500.mat", cellXY)
%
% T       struct returned by ccf_register_fov, or a path to the saved .mat
% cellXY  Nx2 [x y] cell centroids in FOV pixel coordinates, in the same frame the
%         transform was fitted in. If the transform was fitted with flipLR (mirror2),
%         pass "flipApplied", false and this function will flip the x coordinates for you.
%
% Name/value
%   flipApplied  true (default) if cellXY are already in flipped coordinates.
%                Set false to have the flip applied here using T.fov_size.
%   mapFile      override the atlas map path stored in T.
%
% Returns a table: cell index, atlas row/col, AP and ML in mm relative to bregma
% (AP positive = anterior, ML positive = right), hemisphere, area acronym and full name.
%
% Getting cellXY from an EXTRACT output, matching what the analysis scripts already do:
%
%   [Sn, H] = extract_load_ext_spatial(frrFile, nNeurons, 2304);  % sparse pixels x cells
%   nPix   = size(Sn, 1);
%   rowS   = rem((0:nPix-1)', H) + 1;
%   colS   = floor((0:nPix-1)' / H) + 1;
%   wsum   = full(sum(Sn, 1))';
%   cx     = full((rowS' * Sn))' ./ max(wsum, eps);   % this is the "row" coordinate
%   cy     = full((colS' * Sn))' ./ max(wsum, eps);
%   cellXY = [cx, cy];   % check the orientation against your summary image before trusting it

p = inputParser;
p.addParameter("flipApplied", true, @(x) islogical(x) || isnumeric(x));
p.addParameter("mapFile", "", @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

if ischar(T) || isstring(T)
    T = load(char(T));
end
if ~isstruct(T) || ~isfield(T, "tform")
    error("ccf_apply_to_cells:BadTransform", "T must be a ccf_register_fov struct or a path to one.");
end

mapFile = string(opt.mapFile);
if mapFile == ""
    mapFile = string(T.map_file);
end
if exist(mapFile, "file") ~= 2
    % transform may have moved machines; fall back to the copy next to this file
    alt = fullfile(fileparts(mfilename("fullpath")), "allen_dorsal_map.mat");
    if exist(alt, "file") == 2
        mapFile = string(alt);
    else
        error("ccf_apply_to_cells:NoMap", "Atlas map not found: %s", mapFile);
    end
end
A = load(mapFile);

xy = double(cellXY);
if size(xy, 2) ~= 2
    error("ccf_apply_to_cells:BadXY", "cellXY must be Nx2 [x y].");
end

% apply the mirror flip if the caller has not already
if T.flipLR && ~opt.flipApplied
    width = T.fov_size(2);
    xy(:, 1) = width - xy(:, 1) + 1;
end

% FOV pixels -> atlas pixels
ap = transformPointsForward(T.tform, xy);
col = ap(:, 1);
row = ap(:, 2);

% look up the label, guarding the edges
[nR, nC] = size(A.label_img);
ri = round(row);
ci = round(col);
inside = ri >= 1 & ri <= nR & ci >= 1 & ci <= nC;
ids = zeros(size(ri));
lin = sub2ind([nR, nC], ri(inside), ci(inside));
ids(inside) = double(A.label_img(lin));

% id -> acronym / name
acr = strings(size(ids));
nam = strings(size(ids));
acr(:) = "(outside)";
nam(:) = "(outside brain or FOV)";
atlasIds = double(A.ids);
atlasAcr = string(A.acronyms(:));
atlasNam = string(A.names(:));
[tf, loc] = ismember(ids, atlasIds);
acr(tf) = atlasAcr(loc(tf));
nam(tf) = atlasNam(loc(tf));

% stereotaxic-style coordinates relative to bregma
resUm = double(A.res_um);
bregma = double(A.bregma_rc) + 1;      % saved 0-based
ap_mm = -(row - bregma(1)) * resUm / 1000;   % anterior positive
ml_mm =  (col - bregma(2)) * resUm / 1000;   % right positive
hemi = strings(size(ml_mm));
hemi(ml_mm >= 0) = "R";
hemi(ml_mm < 0) = "L";

tbl = table((1:numel(ids))', row, col, ap_mm, ml_mm, hemi, acr, nam, ...
    'VariableNames', {'cell', 'atlas_row', 'atlas_col', 'AP_mm', 'ML_mm', 'hemisphere', ...
                      'area_acronym', 'area_name'});

nOut = sum(~tf);
fprintf("labelled %d cells across %d areas", height(tbl), numel(unique(acr(tf))));
if nOut > 0
    fprintf("  (%d fell outside the atlas)", nOut);
end
fprintf("\n");

% quick summary, most populated first
[u, ~, g] = unique(acr);
cnt = accumarray(g, 1);
[cnt, ord] = sort(cnt, "descend");
u = u(ord);
fprintf("  top areas: ");
for k = 1:min(8, numel(u))
    fprintf("%s %d (%.0f%%)  ", u(k), cnt(k), 100 * cnt(k) / height(tbl));
end
fprintf("\n");
end
