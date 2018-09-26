function [griddedPdf, varargout] = geodata_e_pdf(sData, var, nBins, varargin)

strType = 'separate';

strPtsUse = 'all';
indIn = [];
rng = nan(2,1);
pdfTyp = 'frac';
noVal = 'zero';
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
   elseif regexpbl(varargin{ii}, {'rng', 'range'})
       rng = varargin{ii+1};
   elseif regexpbl(varargin{ii}, 'type')
       pdfTyp = varargin{ii+1};
   elseif regexpbl(varargin{ii}, {'no', 'value'})
       noVal = varargin{ii+1};
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
        %If seperate, range is calculated for each member seperately (if
        %not input)
        
        %Initialize output arrays
        griddedPdf = cell(nMem, 1);
            [griddedPdf{:}] = deal(nan(nLat, nLon, nBins));
        griddedVal = griddedPdf;
        griddedRng = cell(nMem, 1);
            [griddedRng{:}] = deal(nan(nLat, nLon, 2));

        for kk = 1 : nMem
            for ll = 1 : nPts
                [...
                    griddedPdf{kk}(rUse(ll), cUse(ll), :), ...
                    griddedVal{kk}(rUse(ll), cUse(ll), :), ...
                    griddedRng{kk}(rUse(ll), cUse(ll), :)] ...
                    = e_pdf(squeeze(sData{kk}.(var)(:, rUse(ll), cUse(ll))), nBins, 'rng', rng, 'type', pdfTyp, 'no_value', noVal);
            end
        end
	elseif regexpbl(strType, 'com')
        %Calculate range:
        if all(isnan(rng))
            rng = nan([2, [nLat, nLon]], 'single');
            for ll = 1 : nPts
                tempRng = [];
                for kk = 1 : nMem
                    tempRng = [tempRng; squeeze(sData{kk}.(var)(:, rUse(ll), cUse(ll)))];
                end
                rng(:, rUse(ll), cUse(ll)) = [nanmin(tempRng), nanmax(tempRng)];
            end
        end
        szRng = size(rng);
        
        griddedPdf = nan(nLat, nLon, nBins);
        griddedVal = griddedPdf;
        griddedRng = nan(nLat, nLon, 2);
        
        for ll = 1 : nPts
            tempData = [];
            for kk = 1 : nMem
                tempData = [tempData; squeeze(sData{kk}.(var)(:, rUse(ll), cUse(ll)))];
            end
            if isequal(szRng(end-1:end), [nLat, nLon])
                rngCurr = rng(:, rUse(ll), cUse(ll));
            elseif numel(rng) == 2
                rngCurr = rng;
            else
                error('geodataEpdf:rngCurrSize','The range current variable has an unexpected size.');
            end
            [...
                griddedPdf(rUse(ll), cUse(ll), :), ...
                griddedVal(rUse(ll), cUse(ll), :), ...
                griddedRng(rUse(ll), cUse(ll), :)] ...
                = e_pdf(tempData, nBins, 'rng', rngCurr, 'type', pdfTyp, 'no_value', noVal);
        end
    elseif regexpbl(strType, {'avg','average'})
        %Calculate range:
        if all(isnan(rng))
            rng = nan([2, [nLat, nLon]], 'single');
            for ll = 1 : nPts
                tempRng = [];
                for kk = 1 : nMem
                    tempRng = [tempRng; squeeze(sData{kk}.(var)(:, rUse(ll), cUse(ll)))];
                end
                rng(:, rUse(ll), cUse(ll)) = [nanmin(tempRng), nanmax(tempRng)];
            end
        end
        szRng = size(rng);
        
        griddedPdf = nan(nLat, nLon, nBins);
        griddedVal = griddedPdf;
        griddedRng = nan(nLat, nLon, 2);

        for ll = 1 : nPts
            if isequal(szRng(end-1:end), [nLat, nLon])
                rngCurr = rng(:, rUse(ll), cUse(ll));
            elseif numel(rng) == 2
                rngCurr = rng;
            else
                error('geodataEpdf:rngCurrSize','The range current variable has an unexpected size.');
            end
            
            tmpPdf = nan(nMem, nBins);
            tmpVal = nan(nMem, nBins);
            tmpRng = nan(nMem, 2);
            
            for kk = 1 : nMem
                [tmpPdf(kk,:), tmpVal(kk,:), tmpRng(kk,:)] ...
                    = e_pdf(squeeze(sData{kk}.(var)(:, rUse(ll), cUse(ll))), nBins, 'rng', rngCurr, 'type', pdfTyp, 'no_value', noVal);
            end

            griddedPdf(rUse(ll), cUse(ll), :) = nanmean(tmpPdf,1);
            griddedVal(rUse(ll), cUse(ll), :) = nanmean(tmpVal,1);
            griddedRng(rUse(ll), cUse(ll), :) = [nanmin(tmpRng(:)), nanmax(tmpRng(:))];
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
    griddedPdf = nan(nLat, nLon, nBins);
    griddedVal = griddedPdf;
    griddedRng = nan(nLat, nLon, 2);
   
    for ll = 1 : nPts
        [...
            griddedPdf(rUse(ll), cUse(ll), :), ...
            griddedVal(rUse(ll), cUse(ll), :), ...
            griddedRng(rUse(ll), cUse(ll), :)] ...
            = e_pdf(squeeze(sData.(var)(:, rUse(ll), cUse(ll))), nBins, 'rng', rng, 'type', pdfTyp, 'no_value', noVal);
    end
else
    error('geodataECdf:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.']);
end


%Write range as variable output argument
if nargout > 1
   varargout{1} = griddedVal;
   if nargout > 2
        varargout{2} = griddedRng;
   end
end

