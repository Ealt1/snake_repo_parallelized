function [gcell_idx, nvdim] = expandProgressGroup(gpart, ...
    cell_idx,  fx, domain_cell_idx, nvdim, bad_cell_idx)
% Input:
% Output: 
% Assume: 
%   cell_idx is transient under fx


if (nargin < 6)
    bad_cell_idx = [];
end
if (nargin < 5)
    nvdim = setdiff(1:1:gpart.dim, gpart.ddim);
    nvdim = [nvdim, -nvdim];
end
if (nargin < 4)
    domain_cell_idx = gpart.cell_list;
end

% specific tp FCTM
Igrid = gpart.grid{gpart.ddim};
cell_cord = gpart.idx2cord(cell_idx);
if (any(Igrid(cell_cord(gpart.ddim,:))~=fx.u(1)))
    gcell_idx = [];
    nvdim = [];
    return;
end

neighbor_cell_idx = gpart.getNeighbor_idx(cell_idx, 0);
new_cell_idx = setdiff(neighbor_cell_idx, cell_idx);
gcell_idx = cell_idx;
for idx = setdiff(new_cell_idx, bad_cell_idx)
    [~, nvdim] = gpart.isTransient(gcell_idx, fx, nvdim);
    [ts, nvdim] = gpart.isTransient(union(gcell_idx, idx), fx, nvdim);
    if ts
        gcell_idx = union(gcell_idx, idx);
        intersect(gcell_idx, domain_cell_idx);
    else
        bad_cell_idx = union(bad_cell_idx, idx);
    end
end

while (~isempty(setdiff(gcell_idx,cell_idx)) && length(gcell_idx)<=1)
    % gpart.idx2cord(cell_idx)
    % length(gcell_idx)
    % nvdim
    % gpart.idx2cord(gcell_idx)
    cell_idx = gcell_idx;
    [gcell_idx, nvdim] = expandProgressGroup...
        (gpart, cell_idx, fx, domain_cell_idx, nvdim, bad_cell_idx);
end




