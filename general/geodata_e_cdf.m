function griddedCdf = geodata_e_cdf(sData, var, nBins, varargin)

strType = 'separate';

strPtsUse = 'all';
indIn = [];
for ii = 1 : numel(varargin(:))
   if regexpbl(varargin{ii}, {'pts', 'use'}, 'and')
       if isnumeric(varargin{ii+1})
            strPtsUse = 'subset';
            indIn = varargin{ii+1};
       else
            warning('geodataECdf:wrongargType', ...
                ['Optional argument pts_use selected but argument ' ...
                'containing points is type ' class(varargin{ii+1}) ...
                ', which is not expected (must be numeric).'])
       end
   elseif regexpbl(varargin{ii}, 'method')
       strType = varargin{ii+1};
   end
end



if iscell(sData)
    nMem = numel(sData(:));
    
    %Find data size:
    szData = size(sData{1}.(var));
    if numel(szData) == 2
        nLat = szData(1);
        nLon = szData(2);
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
    
    if regexpbl(strType, 'sep')
        %Initialize output arrays
        griddedCdf = cell(nMem, 1);
        [griddedCdf{:}] = deal(nan(nLat, nLon, nBins+1));

        for kk = 1 : nMem
            for ll = 1 : nPts
                [~, griddedCdf{kk}(rUse(ll), cUse(ll), :)] = e_cdf(squeeze(sData{kk}.(var)(:, rUse(ll), cUse(ll))), 'bins', nBins);
            end
        end
	elseif regexpbl(strType, 'com')
        griddedCdf = nan(nLat, nLon, nBins+1);
        for ll = 1 : nPts
            tempData = [];
            for kk = 1 : nMem
                tempData = [tempData; squeeze(sData{kk}.(var)(:, rUse(ll), cUse(ll)))];
            end
            [~, griddedCdf(rUse(ll), cUse(ll), :)] = e_cdf(tempData, 'bins', nBins);
        end
    else
        error('geodataECdf:unknownMethod',['Analysis method ' strType ' has not been programmed for.']);
    end
elseif isstruct(sData)
    %Find data size:
    szData = size(sData.(var));
    if numel(szData) == 2
        nLat = szData(1);
        nLon = szData(2);
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
    
    %Initialize output arrays
    griddedCdf = nan(nLat, nLon, nBins+1);

    for ll = 1 : nPts
        [~, griddedCdf(rUse(ll), cUse(ll), :)] = e_cdf(squeeze(sData.(var)(:, rUse(ll), cUse(ll))), 'bins', nBins);
    end
else
    error('geodataECdf:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.']);
end

