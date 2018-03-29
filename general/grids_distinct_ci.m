function gridCiDiff = grids_distinct_ci(grid1, ci1, grid2, ci2, varargin)

indUse = [];
if ~isempty(varargin)
    for ii = 1 : numel(varargin(:))
        if strcmpi(varargin{ii}, 'ind')
            indUse = varargin{ii+1};
        end
    end
end

if ndims(grid1) ~= 2 || ndims(grid2) ~= 2
   error('gridsDistinctCi:wrongSize','One of the grids does not have 2 dimensions.'); 
end

if isempty(indUse)
   indUse = (1:numel(grid1));
end

[rUse, cUse] = ind2sub(size(grid1), indUse);

%Note: I'm not sure this is appropriate way to calculate confidence
%that means are different:
gridCiDiff = zeros(size(grid1));
for kk = 1 : numel(indUse)
    if grid1(rUse(kk),cUse(kk)) > grid2(rUse(kk),cUse(kk))
        gridCiDiff(rUse(kk),cUse(kk)) = ci1(1, rUse(kk),cUse(kk)) - ci2(2, rUse(kk),cUse(kk));
    else
        gridCiDiff(rUse(kk),cUse(kk)) = ci2(1, rUse(kk),cUse(kk)) - ci1(2, rUse(kk),cUse(kk));
    end
end
clear kk

gridCiDiff(gridCiDiff < 0) = 0;
gridCiDiff(isnan(gridCiDiff)) = 0;

