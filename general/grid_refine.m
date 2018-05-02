function [lonOut, latOut] = grid_refine(lonIn, latIn, nScl, type)
%Ensures that edges of input and output grids line up (centers may not
%align)

if regexpbl(type, 'edge')
    edgeLat = latIn;
    edgeLon = lonIn;
elseif regexpbl(type, 'center')
    edgeLat = box_edg(latIn);
    edgeLon = box_edg(lonIn);
else
    error('gridRefine:unknownType', ['The type ' type ' has not been programmed for.']);
end

edgeLonOut = linspace(min(edgeLon), max(edgeLon), nScl*(numel(edgeLon)-1)+1);
lonOut = edgeLonOut(1:end-1) + 0.5*diff(edgeLonOut);

edgeLatOut = linspace(max(edgeLat), min(edgeLat), nScl*(numel(edgeLat)-1)+1);
latOut = edgeLatOut(1:end-1) + 0.5*diff(edgeLatOut);