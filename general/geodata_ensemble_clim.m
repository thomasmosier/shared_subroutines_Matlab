function ltMean = geodata_ensemble_clim(sData, var)    

if iscell(sData)
    nMod = numel(sData(:));
    
    szData = size(sData{1}.(var));
    
    tempAvg = nan([nMod,  szData(2:3)], 'single');
    nYrsStck = zeros([nMod,  szData(2:3)], 'single');
    for kk = 1 : nMod
        for ll = 1 : szData(1)
            nYrsStck(kk,:,:) = nYrsStck(kk,:,:) + ~isnan(sData{kk}.(var)(ll,:,:));
        end
        clear ll
        tempAvg(kk,:,:) = squeeze(nYrsStck(kk,:,:)).*squeeze(nanmean(sData{kk}.(var), 1));
    end
    clear kk
    ltMean = squeeze(nansum(tempAvg, 1)) ./ squeeze(sum(nYrsStck, 1));
elseif isstruct(sData)
    ltMean = squeeze(nanmean(sData.(var), 1));
else
    error('geodataEnsembleClim:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.']);
end