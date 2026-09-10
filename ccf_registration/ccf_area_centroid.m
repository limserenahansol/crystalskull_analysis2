function out = ccf_area_centroid(acronyms, varargin)
%CCF_AREA_CENTROID  Atlas pixel coordinates of named Allen areas, for functional registration.
%
%   ccf_area_centroid()                       % list every surface area in the dorsal map
%   P = ccf_area_centroid(["VISp","SSp-bfd","SSp-ul"])          % both hemispheres
%   P = ccf_area_centroid(["VISp","SSp-bfd"], "hemisphere", "R")
%
% Returns a table with the area, hemisphere, atlas [x y] pixel coordinates and the AP/ML
% position in mm relative to bregma.
%
% Use it to build the ccfPoints argument of ccf_register_fov in "points" mode. If you run a
% quick sensory-mapping session - visual flash, whisker deflection, forelimb tap - the peak
% of each evoked map in your FOV pairs with the centroid of the matching area here, and that
% gives a functional registration that beats skull landmarks.
%
% Caution: a centroid is not the same thing as a response peak. VISp responses to a full-field
% flash sit near the VISp centroid, but a whisker-evoked peak lands in the barrel column for
% the whisker you deflected, not at the centre of SSp-bfd. Treat these as starting values and
% check the QC overlay.

p = inputParser;
p.addParameter("hemisphere", "both", @(s) ismember(string(s), ["both", "L", "R"]));
p.addParameter("mapFile", fullfile(fileparts(mfilename("fullpath")), "allen_dorsal_map.mat"), ...
    @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

A = load(opt.mapFile);
lab = double(A.label_img);
resUm = double(A.res_um);
bregma = double(A.bregma_rc) + 1;
atlasIds = double(A.ids);
atlasAcr = string(A.acronyms(:));
atlasNam = string(A.names(:));
[~, nC] = size(lab);
midCol = bregma(2);

if nargin < 1 || isempty(acronyms)
    % just print what is available and return the full list
    fprintf("%d areas on the dorsal surface:\n", numel(atlasAcr));
    for k = 1:numel(atlasAcr)
        n = sum(lab(:) == atlasIds(k));
        fprintf("  %-10s %-52s %6d px  %.2f mm2\n", atlasAcr(k), atlasNam(k), n, ...
            n * (resUm / 1000)^2);
    end
    out = table(atlasAcr, atlasNam, 'VariableNames', {'acronym', 'name'});
    return
end

acronyms = string(acronyms(:));
hemis = ["L", "R"];
if string(opt.hemisphere) ~= "both"
    hemis = string(opt.hemisphere);
end

rows = {};
for a = acronyms'
    k = find(atlasAcr == a, 1);
    if isempty(k)
        warning("ccf_area_centroid:Unknown", ...
            "'%s' is not on the dorsal surface. Call with no arguments to list what is.", a);
        continue
    end
    for h = hemis
        m = lab == atlasIds(k);
        if h == "L"
            m(:, ceil(midCol):end) = false;
        else
            m(:, 1:floor(midCol)) = false;
        end
        if ~any(m(:))
            continue
        end
        [rr, cc] = find(m);
        % medoid, so the point is guaranteed to lie inside a concave area
        d = (rr - mean(rr)).^2 + (cc - mean(cc)).^2;
        [~, i] = min(d);
        r = rr(i);
        c = cc(i);
        rows(end + 1, :) = {a, atlasNam(k), h, c, r, ...
            -(r - bregma(1)) * resUm / 1000, (c - bregma(2)) * resUm / 1000, ...
            sum(m(:)) * (resUm / 1000)^2}; %#ok<AGROW>
    end
end

if isempty(rows)
    out = table();
    return
end
out = cell2table(rows, 'VariableNames', ...
    {'acronym', 'name', 'hemisphere', 'x', 'y', 'AP_mm', 'ML_mm', 'area_mm2'});
end
