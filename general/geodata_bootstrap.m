function [avgData, ciData] = geodata_bootstrap(sData, varUse, type, BlkBnds, nBt, sig, varargin)

indUse = [];
if ~isempty(varargin)
    for ii = 1 : numel(varargin(:))
        if strcmpi(varargin{ii}, 'ind')
            indUse = varargin{ii+1};
        end
    end
end

%Check that size of all input arrays is same:
for ii = 1 : numel(sData(:))
    szCurr = size(sData{ii}.(varUse));
    
    if numel(szCurr) ~= 3
        error('geodataBootstrap:not3Dim', ['The current array has ' numel(szCurr) ' dimensions. 3 dimensions are required.'])
    end
    
    if ii == 1
        szData = szCurr;
    else
        if ~isequal(szData(2:3), szCurr(2:3))
           error('geodataBootstrap:diffSpatialDim', 'Not all ensemble members have the same spatial dimensions.');
        end
    end
end
clear ii

if isempty(indUse)
   indUse = (1:prod(szData(2:3)));
end

nBlkSz = randi([min(BlkBnds), max(BlkBnds)], [nBt, 1]); %Number of samples used (ranges from 2 to n-1)

%Calculate weighting factors for ensemble member selection:
wgtData = geodata_ensemble_wgt(sData, varUse);

%Randomly choose which ensemble member to use during each bootstrap loop:
%Weight based on number of years present in each ensemble member
indData = rand_weighted((1:numel(wgtData(:))), wgtData/sum(wgtData), nBt);

%Initialize outputs:
avgData = nan(szData(2:3), 'single');
ciData = nan([2, szData(2:3)], 'single');

%Initialize arrays to use in parfor loop
avgDataTemp = nan([numel(indUse), 1]);
ciDataTemp = nan([2, numel(indUse)], 'single');

valLb = sig; %Lower bound significance
valUb = 1 - sig; %Upper bound significance

if isempty(gcp('nocreate'))
    localCluster = parcluster('local'); %Check number of possible "workers"
    parpool(localCluster.NumWorkers); %Dedicate all available cores to parallel calibration
    
    blClose = 1;
else
    blClose = 0;
end
            
parfor ll = 1 : numel(indUse)
    rCurr = rUse(ll);
    cCurr = cUse(ll);

    btDataTemp = nan([nBt, 1]);

    tsDataTemp = btDataTemp;
    
    %Loop over bootstrap iterations:
    for mm = 1 : nBt
        %Keep only bootstrapped mean
        if strcmpi(type, 'stationary')
            btDataTemp(mm) = mean(stationaryBB(sData{indData(mm)}.(varUse)(:,rCurr,cCurr), 1, nBlkSz(mm)));
        else
            error('geodataBootstrap:unknownType',['Type ' type ' has not been programmed for.']);
        end
    end
    
    %Sort long-term mean calculated from bootstrap iterations
    if strcmpi(varUse, varYrDetr)
        tsDataTemp = sort(ltEnsAvgData(rCurr,cCurr) + btDataTemp);
    elseif strcmpi(varUse, varYr)
        tsDataTemp = sort(btDataTemp);
    else
        error('GissAnalysis:uknownAnalysisVar',['The analysis variable, ' varUse ', does not match the programmed options.']);
    end
    %Remove any nan values:
    tsDataTemp(isnan(tsDataTemp)) = [];

    %Calculate mean from bootstrap:
    avgDataTemp(ll) = squeeze(nanmean(tsDataTemp));

    %Calculate confidendence interval for each long-term mean based on bootstrapped means
    ciDataTemp(:,ll) = tsDataTemp(round([valLb,valUb]*numel(tsDataTemp)));
end
clear ll

if blClose == 1
    %Close dedicated workers used for parallel processing:
    delete(gcp('nocreate'));
end

%Format outputs:
avgData(indUse) = avgDataTemp;
for ll = 1 : numel(indUse)
    ciData(:,rUse(ll), cUse(ll)) = ciDataTemp(:,ll);
end