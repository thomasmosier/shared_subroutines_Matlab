function [avgData, ciData, varargout] = geodata_bootstrap(sData, varUse, type, blkBnds, nBt, sig, varargin)

varDate = 'date';
indUse = [];
nRunning = [];
parClose = 1;
if ~isempty(varargin)
    for ii = 1 : numel(varargin(:))
        if strcmpi(varargin{ii}, 'ind')
            indUse = varargin{ii+1};
        elseif regexpbl(varargin{ii}, {'run', 'mean'}, 'and')
            nRunning = varargin{ii+1};
        elseif regexpbl(varargin{ii}, {'par', 'close'}, 'and')
            parClose = varargin{ii+1};
        end
    end
end

if isstruct(sData)
   sData = {sData}; 
end

%Check that size of all input arrays is same:
for ii = 1 : numel(sData(:))
    szCurr = size(sData{ii}.(varUse));
    
    if ii == 1
        szData = szCurr;
    else
        if ~isequal(szData(2:end), szCurr(2:end))
           error('geodataBootstrap:diffSpatialDim', 'Not all ensemble members have the same spatial dimensions.');
        end
    end
end
clear ii

if isempty(indUse)
   indUse = (1:prod(szData(2:end)));
end


nBlkSz = randi([min(blkBnds), max(blkBnds)], [nBt, 1]); %Number of samples used (ranges from 2 to n-1)

%Calculate weighting factors for ensemble member selection:
wgtData = geodata_ensemble_wgt(sData, varUse);

%Randomly choose which ensemble member to use during each bootstrap loop:
%Weight based on number of years present in each ensemble member
indData = rand_weighted((1:numel(wgtData(:))), wgtData/sum(wgtData), nBt);

nMem = numel(sData(:));

valLb = sig; %Lower bound significance
valUb = 1 - sig; %Upper bound significance


%Define all variables that occur in parallel loop (even if not used with specific options set)  
datesLp = [];
runAvgTemp = [];
ciRunningTemp = [];
avgDataTemp = [];
ciDataTemp = [];
if ~isempty(nRunning)
    datesLp = sData{1}.(varDate);
    if numel(sData(:)) > 1
        for ii = 2 : nMem
            datesLp = union(datesLp, sData{ii}.(varDate), 'rows');
        end
    end
    
    nTime = numel(datesLp(:,1));
    nStep = ceil(nRunning/2 - 1);
    
    %Initialize outputs:
    avgData = nan(szData, 'single');
    ciData = nan([2, szData], 'single');

    %Initialize arrays to use in parfor loop
    runAvgTemp = nan([nTime, numel(indUse)]);
    ciRunningTemp = nan([2, nTime, numel(indUse)], 'single');
else
    nTime = 0;
    nStep = nan;
    
    %Initialize outputs:
    avgData = nan(szData, 'single');
    ciData = nan([2, szData], 'single');

    %Initialize arrays to use in parfor loop
    avgDataTemp = nan([numel(indUse), 1]);
    ciDataTemp = nan([2, numel(indUse)], 'single');
end


nDim = numel(szData);
if nDim == 3
    [rUse, cUse] = ind2sub(szData(2:3), indUse);
else
    rUse = szData(2);
    cUse = rUse;
end


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

    if ~isempty(nRunning)
        dataTemp = nan([nTime, nMem], 'single');
        for mm = 1 : nMem
            indT = find(ismember(sData{mm}.(varDate), datesLp, 'rows'));
            if nDim == 3
                dataTemp(indT(1):indT(end),mm) = sData{mm}.(varUse)(:,rCurr,cCurr);
            else
                dataTemp(indT(1):indT(end),mm) = sData{mm}.(varUse)(:,rCurr);
            end
        end
        
        %Calculate running mean for current point (average between all
        %ensemble members
        runAvgCurr = nan([nTime, nMem], 'single');
        for mm = 1 : nMem
            runAvgCurr(:,mm) = runmean(dataTemp(:,mm), nStep, [], 'edge');
        end
        runAvgTemp(:,ll) = nanmean(runAvgCurr, 2);
        
        %Calculate bootstrapped confidence intervals for each time:    
        for zz = 1 : nTime
            btDataTemp = nan([nBt, 1], 'single');
            
            indT = [max(zz - nStep, 1), min(zz + nStep, nTime)];
            %Loop over bootstrap iterations:
            for mm = 1 : nBt
                %Keep only bootstrapped mean
                if strcmpi(type, 'stationary')
                    btDataTemp(mm) = runAvgCurr(zz,indData(mm)) + mean(stationaryBB(detrend(dataTemp(indT(1):indT(2),indData(mm))), 1, nBlkSz(mm)));
                else
                    error('geodataBootstrap:unknownType',['Type ' type ' has not been programmed for.']);
                end
            end
            
            %Sort bootstrapped values
            btDataTemp = sort(btDataTemp);

            %Remove any nan values:
            btDataTemp(isnan(btDataTemp)) = [];

            %Calculate confidence intervals from bootstrap:
            if ~isempty(btDataTemp)
                %Calculate confidendence interval for each long-term mean based on bootstrapped means
                ciRunningTemp(:,zz,ll) = btDataTemp(round([valLb,valUb]*numel(btDataTemp)));
            else
                ciRunningTemp(:,zz,ll) = nan(2,1);
            end
        end
    else
        btDataTemp = nan([nBt, 1]);
            
        %Loop over bootstrap iterations:
        for mm = 1 : nBt
            %Keep only bootstrapped mean
            if strcmpi(type, 'stationary')
                if nDim == 3
                    btDataTemp(mm) = mean(stationaryBB(sData{indData(mm)}.(varUse)(:,rCurr,cCurr), 1, nBlkSz(mm)));
                else
                    btDataTemp(mm) = mean(stationaryBB(sData{indData(mm)}.(varUse)(:,rCurr), 1, nBlkSz(mm)));
                end
            else
                error('geodataBootstrap:unknownType',['Type ' type ' has not been programmed for.']);
            end
        end
        
        %Sort bootstrapped values
        btDataTemp = sort(btDataTemp);

        %Remove any nan values:
        btDataTemp(isnan(btDataTemp)) = [];

        %Calculate mean from bootstrap:
        if ~isempty(btDataTemp)
            avgDataTemp(ll) = squeeze(nanmean(btDataTemp));

            %Calculate confidendence interval for each long-term mean based on bootstrapped means
            ciDataTemp(:,ll) = btDataTemp(round([valLb,valUb]*numel(btDataTemp)));
        else
            avgDataTemp(ll) = nan;
            ciDataTemp(:,ll) = nan(2,1);
        end
    end
end
clear ll

if blClose == 1 && parClose == 1
    %Close dedicated workers used for parallel processing:
    delete(gcp('nocreate'));
end

%Format outputs:
if ~isempty(nRunning)
    %Assign outputs to grids (must account for time dimension in
    %assignment)
    for ll = 1 : numel(indUse)
        avgData(:,indUse(ll)) = runAvgTemp(:,ll);
    end
    
    if nDim == 3
        for ll = 1 : numel(indUse)
            ciData(:,:,rUse(ll), cUse(ll)) = ciRunningTemp(:,:,ll);
        end
    else
        if numel(indUse) == 1
            ciData = ciRunningTemp;
        else
            for ll = 1 : numel(indUse)
                ciData(:,:,rUse(ll)) = ciRunningTemp(:,:,ll);
            end
        end
    end
else
    avgData(indUse) = avgDataTemp;
    if nDim == 3
        for ll = 1 : numel(indUse)
            ciData(:,rUse(ll), cUse(ll)) = ciDataTemp(:,ll);
        end
    else
        for ll = 1 : numel(indUse)
            ciData(:,rUse(ll)) = ciDataTemp(:,ll);
        end
    end
end