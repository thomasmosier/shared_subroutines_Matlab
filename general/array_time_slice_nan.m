function array = array_time_slice_nan(array, tInd, spIndNan)

szArray = size(array);

if numel(szArray) == 3
    [rNan, cNan] = ind2sub(szArray(2:3), spIndNan);
elseif numel(szArray) == 2
    [rNan, cNan] = ind2sub(szArray(1:2), spIndNan);
end

if numel(tInd) == 1
    tInd = repmat(tInd, numel(rNan), 1);
else
    error('printModelPt:multTimeInd',['There are ' num2str(tInd) ' input time indices. Only one is allowed.']);
end
indNan3d = sub2ind(szArray, tInd, rNan(:), cNan(:));

array(indNan3d) = nan;