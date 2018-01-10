function [matNorth, matEast, matZone] = deg2utm_grid(vecLat, vecLon)

if ((any(size(vecLon) == 1) && numel(vecLon) > 1) && (~any(size(vecLat) == 1) && numel(vecLat) > 1) ) || ((~any(size(vecLon) == 1) && numel(vecLon) > 1) && (any(size(vecLat) == 1) && numel(vecLat) > 1) )
    error('deg2utm_grid:latLonVecDim', 'One coordinate input is a 1D vector and the other is a 2D array. Dimensions of these inputs must be the same.');
end



% matNorth1 = nan([numel(vecLat), numel(vecLon)]);
% matEast1 = matNorth; 
% matZone1 = cell(size(matEast));

%Accomodate inputs as vectors or 2D arrays
if any(size(vecLat)) == 1 && any(size(vecLon)) == 1
    vecLon = vecLon(:)';
    vecLat = vecLat(:);
    
    %Initialize outputs:
    matNorth = nan([numel(vecLat), numel(vecLon)]);
    matEast = matNorth; 
    matZone = cell(size(matEast));

    for ii = 1 : numel(vecLat(:,1))
        for jj = 1 : numel(vecLon(1,:))
            [matEast(ii,jj), matNorth(ii,jj), matZone{ii,jj}] ...
                = deg2utm(vecLat(ii), vecLon(jj)); %[x,y,utmzone] = deg2utm(Lat,Lon)
    %         [matEast1(ii,jj), matNorth1(ii,jj), matZone1{ii,jj}] = ll2utm([vecLat(ii), vecLon(jj)]);
        end
    end
else %2D array
    %Initialize outputs:
    matNorth = nan([numel(vecLat(:,1)), numel(vecLon(1,:))]);
    matEast = matNorth; 
    matZone = cell(size(matEast));
    
    for ii = 1 : numel(vecLat(:,1))
        for jj = 1 : numel(vecLon(1,:))
            [matEast(ii,jj), matNorth(ii,jj), matZone{ii,jj}] ...
                = deg2utm(vecLat(ii,jj), vecLon(ii,jj)); %[x,y,utmzone] = deg2utm(Lat,Lon)
    %         [matEast1(ii,jj), matNorth1(ii,jj), matZone1{ii,jj}] = ll2utm([vecLat(ii), vecLon(jj)]);
        end
    end
end


%[east, north, zone] = deg2utm(vecLat(1), vecLon(2))