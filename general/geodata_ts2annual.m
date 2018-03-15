function sData = geodata_ts2annual(sData, varOut, varYrData, type)

varDate = 'date';
varYrDate = [varDate, varYrData];

if iscell(sData)
    for kk = 1 : numel(sData(:))
        sData{kk}.(varYrDate) = unique(sData{kk}.(varDate)(:,1));
        
        szData = size(sData{kk}.(varOut));

        sData{kk}.(varYrData) = nan([numel(sData{kk}.(varYrDate)), szData(2:3)], 'single');
        for ll = 1 : numel(sData{kk}.(varYrDate))
            indCurr = find(sData{kk}.(varDate)(:,1) == sData{kk}.(varYrDate)(ll));

            if isempty(indCurr) 
                sData{kk}.(varYrData)(ll,:,:) = nan;
            else
                if regexpbl(type, {'mean','avg'})
                    sData{kk}.(varYrData)(ll,:,:) = squeeze(nanmean(sData{kk}.(varOut)(indCurr,:,:), 1));
                elseif regexpbl(type, 'sum')
                    sData{kk}.(varYrData)(ll,:,:) = squeeze(nansum(sData{kk}.(varOut)(indCurr,:,:), 1));
                else
                    error('geodataTs2Annual:unknownAggType',['The aggregation type ' type ' has not been programmed for.']);
                end
            end
        end
        clear ll
    end
    clear kk
elseif isstruct(sData)
    sData.(varYrDate) = unique(sData.(varDate)(:,1));
    
    szData = size(sData.(varOut));

    sData.(varYrData) = nan([numel(sData.(varYrDate)), szData(2:3)], 'single');
    for ll = 1 : numel(sData.(varYrDate))
        indCurr = find(sData.(varDate)(:,1) == sData.(varYrDate)(ll));

        if isempty(indCurr) 
            sData.(varYrData)(ll,:,:) = nan;
        else
            if regexpbl(type, {'mean','avg'})
                sData.(varYrData)(ll,:,:) = squeeze(nanmean(sData.(varOut)(indCurr,:,:), 1));
            elseif regexpbl(type, 'sum')
                sData.(varYrData)(ll,:,:) = squeeze(nansum(sData.(varOut)(indCurr,:,:), 1));
            else
                error('geodataTs2Annual:unknownAggType',['The aggregation type ' type ' has not been programmed for.']);
            end
        end
    end
    clear ll
else
    error('geodataOceanNan:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.']);
end