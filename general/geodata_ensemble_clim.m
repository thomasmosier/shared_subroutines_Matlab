function varargout = geodata_ensemble_clim(sData, var, varargin)
indUse = [];
blStructOut = 0;
if ~isempty(varargin)
    for ii = 1 : numel(varargin(:))
        if strcmpi(varargin{ii}, 'ind')
            indUse = varargin{ii+1};
        elseif regexpbl(varargin{ii}, {'type', 'each'}, 'and')
            blStructOut = 1;
        end
    end
end

varOut = [var '_avg'];

if iscell(sData)
    nMod = numel(sData(:));
    
    szData = size(sData{1}.(var));
    
    if blStructOut == 0
        if numel(szData) == 3
            tempAvg = nan([nMod,  szData(2:3)], 'single');
            nYrsStck = zeros([nMod,  szData(2:3)], 'single');
            for kk = 1 : nMod
                nTimeCurr = numel(sData{kk}.(var)(:,1,1));
                for ll = 1 : nTimeCurr
                    nYrsStck(kk,:,:) = nYrsStck(kk,:,:) + ~isnan(sData{kk}.(var)(ll,:,:));
                end
                clear ll
                tempAvg(kk,:,:) = squeeze(nYrsStck(kk,:,:)).*squeeze(nanmean(sData{kk}.(var), 1));
            end
            clear kk
            varargout{1} = squeeze(nansum(tempAvg, 1)) ./ squeeze(sum(nYrsStck, 1));

            if ~isempty(indUse)
                varargout{1}(setdiff((1:numel(varargout{1})), indUse)) = nan;
            end
        elseif numel(szData) == 2
            tempAvg = nan([nMod,  szData(2)], 'single');
            nYrsStck = zeros([nMod,  szData(2)], 'single');
            for kk = 1 : nMod
                nTimeCurr = numel(sData{kk}.(var)(:,1));
                for ll = 1 : nTimeCurr
                    nYrsStck(kk,:,:) = nYrsStck(kk,:,:) + ~isnan(sData{kk}.(var)(ll,:));
                end
                clear ll
                tempAvg(kk,:,:) = squeeze(nYrsStck(kk,:,:)).*squeeze(nanmean(sData{kk}.(var), 1));
            end
            clear kk
            varargout{1} = squeeze(nansum(tempAvg, 1)) ./ squeeze(sum(nYrsStck, 1));
        else
            error('geodataEnsembleClim:dimWrong','The input data does not have an appropriate number of dimensions.');
        end
    elseif blStructOut == 1
        for kk = 1 : nMod
            sData{kk}.(varOut) = squeeze(nanmean(sData{kk}.(var), 1));
            
            if numel(szData) == 3
                sData{kk}.(varOut)(setdiff((1:numel(sData{kk}.(varOut))), indUse)) = nan;
            end
        end
        clear kk
        
        varargout{1} = sData;
    end
elseif isstruct(sData)
    if blStructOut == 0
        varargout{1} = squeeze(nanmean(sData.(var), 1));
        
        if ~isempty(indUse)
            varargout{1}(setdiff((1:numel(varargout{1})), indUse)) = nan;
        end
    elseif blStructOut == 1
        sData.(varOut) = squeeze(nanmean(sData.(var), 1));
        
        sData.(varOut)(setdiff((1:numel(sData.(varOut))), indUse)) = nan;
        
        varargout{1} = sData;
    end
else
    error('geodataEnsembleClim:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.']);
end