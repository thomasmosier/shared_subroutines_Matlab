function sData = geodata_anomaly(sData, varUse, varOut, nAvg, type, varargin)

strPtsUse = 'all';
indIn = [];
for ii = 1 : numel(varargin(:))
   if regexpbl(varargin{ii}, {'pts', 'use'}, 'and')
       if isnumeric(varargin{ii+1})
            strPtsUse = 'subset';
            indIn = varargin{ii+1};
       else
            warning('geodataAnomaly:wrongargType', ...
                ['Optional argument pts_use selected but argument ' ...
                'containing points is type ' class(varargin{ii+1}) ...
                ', which is not expected (must be numeric).'])
       end
   end
end

if iscell(sData)
    nMem = numel(sData(:));
    
    szData = size(sData{1}.(varUse));
    if numel(szData) == 2
        error('geodataAnomaly:twoD','The input grid is 2d but needs to be 3d.');
    elseif numel(szData) == 3
        nLat = szData(2);
        nLon = szData(3);
    end
    
    %Find indices to use
    if strcmpi(strPtsUse, 'all')
        indUse = (1:nLat*nLon);
    elseif strcmpi(strPtsUse, 'subset')
        indUse = indIn;
    end
    nPts = numel(indUse);
    
    [rUse, cUse] = ind2sub([nLat, nLon], indUse);
    
    for kk = 1 : nMem
        %Initialize output arrays
        sData{kk}.(varOut) = nan(szData, 'single');

        for ll = 1 : nPts
            if regexpbl(type, {'add','sub'})
                sData{kk}.(varOut)(:, rUse(ll), cUse(ll)) = sData{kk}.(varUse)(:, rUse(ll), cUse(ll)) - runmean(sData{kk}.(varUse)(:, rUse(ll), cUse(ll)), ceil(nAvg/2 - 1), [], 'edge');
            elseif regexpbl(type, {'mult','frac'})
                sData{kk}.(varOut)(:, rUse(ll), cUse(ll)) = sData{kk}.(varUse)(:, rUse(ll), cUse(ll)) ./ runmean(sData{kk}.(varUse)(:, rUse(ll), cUse(ll)), ceil(nAvg/2 - 1), [], 'edge');
            else
                error('geodataAnomaly:unknownType',['The type ' type ' has not been programmed for.']);
            end
        end
    end
    
elseif isstruct(sData)
    %Find data size:
    szData = size(sData.(varUse));
    if numel(szData) == 2
        error('geodataAnomaly:twoD','The input grid is 2d but needs to be 3d.');
    elseif numel(szData) == 3
        nLat = szData(2);
        nLon = szData(3);
    end
    
    %Find indices to use
    if strcmpi(strPtsUse, 'all')
        indUse = (1:nLat*nLon);
    elseif strcmpi(strPtsUse, 'subset')
        indUse = indIn;
    end
    nPts = numel(indUse);
    
    [rUse, cUse] = ind2sub([nLat, nLon], indUse);
      
    for ll = 1 : nPts
        if regexpbl(type, {'add','sub'})
            sData.(varOut)(:, rUse(ll), cUse(ll)) = sData.(varUse)(:, rUse(ll), cUse(ll)) - runmean(sData.(varUse)(:, rUse(ll), cUse(ll)), ceil(nAvg/2 - 1), [], 'edge');
        elseif regexpbl(type, {'mult','frac'})
            sData.(varOut)(:, rUse(ll), cUse(ll)) = sData.(varUse)(:, rUse(ll), cUse(ll)) ./ runmean(sData.(varUse)(:, rUse(ll), cUse(ll)), ceil(nAvg/2 - 1), [], 'edge');
        else
            error('geodataAnomaly:unknownType',['The type ' type ' has not been programmed for.']);
        end
    end 
  
else
    error('geodataAnomaly:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.']);
end