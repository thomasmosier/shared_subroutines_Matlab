function sOut = geodata_spatial_avg(sData, varUse, lonUse, latUse)

varLon = 'longitude';
varLat = 'latitude';
varDate = 'date';

if iscell(sData)
    nMem = numel(sData(:));
    
    sOut = cell(nMem, 1);
    for kk = 1 : nMem
        sOut{kk} = struct;
        sOut{kk}.(varDate) = sData{kk}.(varDate);
        
        indLon = [find(sData{kk}.(varLon) >= min(lonUse), 1, 'first'), find(sData{kk}.(varLon) <= max(lonUse), 1, 'last')];
        indLat = [find(sData{kk}.(varLat) <= max(latUse), 1, 'first'), find(sData{kk}.(varLat) >= min(latUse), 1, 'last')];
        
        area = area_geodata(sData{kk}.(varLon)(indLon(1):indLon(end)),sData{kk}.(varLat)(indLat(1):indLat(end)), 'c');
        sumArea = sum2d(area);
        nTs = numel(sData{kk}.(varUse)(:,1,1));
        
        if ~isempty(indLon) && ~isempty(indLat)
            sOut{kk}.(varUse) = nan(nTs, 1);
            
            for jj = 1 : nTs
                sOut{kk}.(varUse)(jj) = sum2d(area.*squeeze(sData{kk}.(varUse)(jj,indLat(1):indLat(end),indLon(1):indLon(end))))/sumArea;
            end
            clear jj
        else
            error('geodataSpatialAvg:outsideBnds', 'The input data does not include the specified region.');
        end
    end
    clear kk
elseif isstruct(sData)
    indLon = [find(sData.(varLon) >= min(lonUse), 1, 'first'), find(sData.(varLon) <= max(lonUse), 1, 'first')];
    indLat = [find(sData.(varLat) <= max(latUse), 1, 'first'), find(sData.(varLat) >= min(latUse), 1, 'first')];

    area = area_geodata(sData.(varLon)(indLon(1):indLon(end)),sData.(varLat)(indLat(1):indLat(end)), 'c');
    sumArea = sum2d(sumArea);
    nTs = numel(sData.(varUse)(:,1,1));

    if ~isempty(indLon) && ~isempty(indLat)
        sOut = struct;
        sOut.(varDate) = sData.(varDate);

        sOut.(varUse) = nan(nTs, 1);
        for jj = 1 : nTs
            sOut.(varUse)(jj) = sum2d(area.*squeeze(sData.(varUse)(jj,indLat(1):indLat(end),indLon(1):indLon(end))))/sumArea;
        end
        clear jj
    else
        error('geodataSpatialAvg:outsideBnds', 'The input data does not include the specified region.');
    end
else
    error('geodataSpatialAvg:unknownClass', ['The input data are of class' class(sData) ', which is not programmed for.']);
end