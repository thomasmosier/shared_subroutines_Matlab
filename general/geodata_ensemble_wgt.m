function wgt = geodata_ensemble_wgt(sData, var)

if iscell(sData)
    wgt = nan([numel(sData), 1]);

    for ii = 1 : numel(wgt)
        if ismatrix(sData{ii}.(var))
            wgt(ii) = numel(sData{ii}.(var)(:,1));
        elseif ndims(sData{ii}.(var)) == 3
            wgt(ii) = numel(sData{ii}.(var)(:,1,1));
        end
    end
    clear ii
elseif isstruct(sData)
    wgt = 1;
else
    error('geodataEnsembleWgt:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.']);
end