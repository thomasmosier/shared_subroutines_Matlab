function [val, ind1, ind2, ind3] = min3d(array)

nInd1 = numel(array(:,1,1));
tempArray = nan([nInd1, 1]);

for ii = 1 : numel(array(:,1,1))
    [tempArray(ii), ~, ~] = min2d(squeeze(array(ii,:,:)));
end

[~, ind1] = min(tempArray);

[val, ind2, ind3] = min2d(array(ind1, :, :));