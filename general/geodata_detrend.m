function sData = geodata_detrend(sData, varIn, varOut)


if iscell(sData)
    for kk = 1 : numel(sData(:))
        szData = size(sData{kk}.(varIn));
        
        sData{kk}.(varOut) = nan(szData, 'single');
        
        indLp = find(~isnan(sData{kk}.(varIn)(1,:,:)));
        [rLp, cLp] = ind2sub(szData(2:end), indLp);

        for ll = 1 : numel(indLp)
            sData{kk}.(varOut)(:,rLp(ll),cLp(ll)) = detrend(squeeze(sData{kk}.(varIn)(:,rLp(ll),cLp(ll))));
        end
        clear ll
    end
    clear kk
elseif isstruct(sData)
        szData = size(sData.(varIn));
        
        sData.(varOut) = nan(szData, 'single');
        
        indLp = find(~isnan(sData.(varIn)(1,:,:)));
        [rLp, cLp] = ind2sub(szData(2:end), indLp);

        for ll = 1 : numel(indLp)
            sData.(varOut)(:,rLp(ll),cLp(ll)) = detrend(squeeze(sData.(varIn)(:,rLp(ll),cLp(ll))));
        end
        clear ll
else
    error('geodataDetrend:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.']);
end