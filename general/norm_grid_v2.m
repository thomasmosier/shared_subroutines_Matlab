function [dataCor, wgtIn, wgtOut] = norm_grid_v2(dataIn, areaIn, latIn, lonIn, dataOut, areaOut, latOut, lonOut, latBnds, lonBnds, strNorm)


%Initialize corrected output
dataCor = nan(size(dataOut));

%Crop if lat/lon bounds provided (NOT NAN)
if all(~isnan(latBnds)) && all(~isnan(lonBnds))
    indLatIn = find(latIn >= min(latBnds) & latIn <= max(latBnds));
    indLonIn = find(lonIn >= min(lonBnds) & lonIn <= max(lonBnds));

    [indLonInMesh, indLatInMesh] = meshgrid(indLonIn, indLatIn);
    indIn = sub2ind(size(dataIn), indLatInMesh(:), indLonInMesh(:));
    
    indLatOut = find(latOut >= min(latBnds) & latOut <= max(latBnds));
    indLonOut = find(lonOut >= min(lonBnds) & lonOut <= max(lonBnds));

    [indLonOutMesh, indLatOutMesh] = meshgrid(indLonOut, indLatOut);
    indOut = sub2ind(size(dataOut), indLatOutMesh(:), indLonOutMesh(:));
    
    dataIn = reshape(dataIn(indIn), [numel(indLatIn), numel(indLonIn)]);
    areaIn = reshape(areaIn(indIn), [numel(indLatIn), numel(indLonIn)]);
    latIn  = latIn(indLatIn);
    lonIn  = lonIn(indLonIn);
    
    dataOut = reshape(dataOut(indOut), [numel(indLatOut), numel(indLonOut)]);
    areaOut = reshape(areaOut(indOut), [numel(indLatOut), numel(indLonOut)]);
    latOut = latOut(indLatOut);
    lonOut = lonOut(indLonOut);
else
    indLatOut = (1:numel(latOut));
    indLonOut = (1:numel(lonOut));
    
    [indLonOutMesh, indLatOutMesh] = meshgrid(indLonOut, indLatOut);
    indOut = sub2ind(size(dataOut), indLatOutMesh(:), indLonOutMesh(:));
end


%Find common lat/lon coordinates:
edgLonIn  = box_edg(lonIn);
edgLonOut = box_edg(lonOut);

edgLatIn  = box_edg(latIn);
edgLatOut = box_edg(latOut);

iLonIn  = nan(1,2);
iLonOut = nan(1,2);
iLatIn  = nan(1,2);
iLatOut = nan(1,2);

if edgLonIn(1) < edgLonOut(1)
    iLonOut(1) = 1;
    [~, iLonIn(1)] = min(abs(edgLonIn - edgLonOut(1)));
else
    iLonIn(1) = 1;
    [~, iLonOut(1)] = min(abs(edgLonOut - edgLonIn(1)));
end

if edgLatIn(1) > edgLatOut(1)
    iLatOut(1) = 1;
    [~, iLatIn(1)] = min(abs(edgLatIn - edgLatOut(1)));
else
    iLatIn(1) = 1;
    [~, iLatOut(1)] = min(abs(edgLatOut - edgLatIn(1)));
end

if edgLonIn(end) > edgLonOut(end)
    iLonOut(2) = numel(edgLonOut)-1;
    [~, indLonInTemp] = min(abs(edgLonIn - edgLonOut(end)));
    iLonIn(2) = indLonInTemp - 1;
else
    iLonIn(2) = numel(edgLonIn)-1;
    [~, indLonOutTemp] = min(abs(edgLonOut - edgLonIn(end)));
    iLonOut(2) = indLonOutTemp - 1;
end

if edgLatIn(end) < edgLatOut(end)
    iLatOut(2) = numel(edgLatOut)-1;
    [~, indLatInTemp] = min(abs(edgLatIn - edgLatOut(end)));
    iLatIn(2) = indLatInTemp - 1;
else
    iLatIn(2) = numel(edgLatIn)-1;
    [~, indLatOutTemp] = min(abs(edgLatOut - edgLatIn(end)));
    iLatOut(2) = indLatOutTemp - 1;
end


wgtIn  = sum2d( areaIn( iLatIn(1): iLatIn(2), iLonIn(1): iLonIn(2)).* dataIn( iLatIn(1): iLatIn(2), iLonIn(1): iLonIn(2)))/sum2d( areaIn( iLatIn(1): iLatIn(2), iLonIn(1): iLonIn(2)));
wgtOut = sum2d(areaOut(iLatOut(1):iLatOut(2),iLonOut(1):iLonOut(2)).*dataOut(iLatOut(1):iLatOut(2),iLonOut(1):iLonOut(2)))/sum2d(areaOut(iLatOut(1):iLatOut(2),iLonOut(1):iLonOut(2)));

if isnan(wgtIn)
    warning('normGrid:wgtInNan','No normalization is being applied because the input weight is nan.');
    return
end

if isnan(wgtOut)
    warning('normGrid:wgtOutNan','No normalization is being applied because the output weight is nan.');
    return
end


if regexpbl(strNorm, {'wgt','mult'}, 'and')
%     norm = (mean2d(dataOut)) / (mean2d(dataIn));
    norm = wgtOut / wgtIn;
    
    if isnan(norm)
       warning('normGrid:normNan', 'The normalization factor is nan. It is therefore being set to 1.');
       norm = 1;
    end
    dataCor(indOut) = dataOut(:) / norm;

    if norm > 10 || norm < 0.1
       warning('normGrid:GT10', ['The normalization factor is ' num2str(norm) ', which is larger than expected.']) 
    end
elseif regexpbl(strNorm, {'wgt','add'}, 'and')
    norm = wgtOut - wgtIn;
    
    dataCor(indOut) = dataOut(:) - norm;

    if norm > 5 || norm < -5
       warning('interp2Norm:GT5', ['The normalization factor is ' num2str(norm) ', which is larger than expected.']) 
    end
    
elseif regexpbl(strNorm, {'none'})
    dataCor(indOut) = dataOut(:);
else
    error('interpNorm:unknownType', ['The normalization method ' ...
        strNorm ' has not been programmed for.']);
end
