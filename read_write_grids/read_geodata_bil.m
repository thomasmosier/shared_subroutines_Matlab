function sData = read_geodata_bil(pathData, varReq, lonLd, latLd, yrsLd, mnthsLd)

varLon = 'longitude';
varLat = 'latitude';
varDate = 'date';


if regexpbl(pathData, 'PRISM')
    if regexpbl(varReq,{'tmp','tas','tmean'})   %Necessary because WorldClim uses 'tmean' instead of 'tmn'.
        varLd = 'tmean';
    elseif regexpbl(varReq,{'tmn','tmin','tasmin'}) 
        varLd = 'tmin';
    elseif regexpbl(varReq,{'tmx','tmax','tasmax'}) 
        varLd = 'tmax';
    elseif regexpbl(varReq,'pre') || strcmpi(varReq,'pr')
        varLd = 'ppt';
    else
        error('readGeodataBil:unknownWCvar',[varReq ' is an '...
            'unknown variable for the PRISM dataset. Program for this if correct.']);
    end
    
    type = 'PRISM';
elseif regexpbl(pathData, {'wc','worldClim'})
    if regexpbl(varReq,{'tmp','tas','tmean'})   %Necessary because WorldClim uses 'tmean' instead of 'tmn'.
        varLd = 'tmean';
    elseif regexpbl(varReq,{'tmn','tmin','tasmin'})
        varLd = 'tmin';
    elseif regexpbl(varReq,{'tmx','tmax','tasmax'})
        varLd = 'tmax';
    elseif regexpbl(varReq,'pre') || strcmpi(varReq,'pr')
        varLd = 'prec';
    elseif regexpbl(varReq,{'alt','elev','DEM'})
        varLd = 'alt';
    else
        error('readGeodataBil:unknownWCvar',[varReq ' is an '...
            'unknown variable for the WorldClim dataset. Program for this if correct.']);
    end
    
    type = 'WorldClim';
elseif regexpbl(pathData, 'APHRO')
    if regexpbl(varReq,{'tmp','tas','tmean'})   %Necessary because WorldClim uses 'tmean' instead of 'tmn'.
        varLd = 'tave';
    elseif regexpbl(varReq,'pre') || strcmpi(varReq,'pr')
        varLd = 'precip';
    else
        error('readGeodataBil:unknownAPHROvar',[varReq ' is an '...
            'unknown variable for the APHRODITE dataset. Program for this if correct.']);
    end
    
    type = 'APHRODITE';
else
    error('readGeodataBil:unknownSource', ['The file ' pathData ' does not match any programmed binary file sources (WorldClim, PRISM, or APHRODITE).'])
end

%Initialize output:
sData = struct;
sData.attTime = cell(0,2);

%Check if WordClim path includes file or is general directory
[root,file,extTemp] = fileparts(pathData);
%If includes file, remove so that it searches directory 
%(necessary because there are both hdr and data files):
if ~isempty(extTemp)
    pathData = root;
else
    pathData = fullfile(root,file);
end

%Check naming and find file for current month:
if strcmpi(type, 'PRISM')
    if mnthsLd < 10
        mnthStr = ['0' num2str(mnthsLd)];
    else
        mnthStr = num2str(mnthsLd);
    end

    filePrismHdr = dir(fullfile(pathData, strcat('PRISM_', varLd, '*', mnthStr, '_bil.hdr')));
    filePrismData = dir(fullfile(pathData, strcat('PRISM_', varLd, '*', mnthStr, '_bil.bil')));
    
    filePrismData = extract_field(filePrismData, 'name');

    %Assign date info based on whether files are time-series or not
    indUnd = regexpi(filePrismData{1}, '_');
    numTest = filePrismData{1}(indUnd(end-1)+1 : indUnd(end)-1);
    if numel(numTest) == 2
        prismType = 'clim';
    elseif numel(numTest) == 6
        prismType = 'ts';

        dateAll = nan(numel(filePrismData),2);
        for zz = 1 : numel(dateAll(:,1))
            indUnd = regexpi(filePrismData{zz}, '_');
            dateAll(zz,1) = str2double(filePrismData{zz}(indUnd(end-1)+1 : indUnd(end-1)+4));
            dateAll(zz,2) = str2double(filePrismData{zz}(indUnd(end-1)+5 : indUnd(end)-1));
        end

        %Select which to keep:

        indYrKeep = find(dateAll(:,1) >= yrsLd(1) & dateAll(:,1) <= yrsLd(2));
        indMnthKeep = find(dateAll(:,2) == mnthsLd);
        indKeep = intersect(indYrKeep, indMnthKeep);
        if numel(indKeep) > 1
            filePrismData = filePrismData(indKeep);
        else
            error('readGeodataBil:noPrism','No PRISM files were found for the requested dates.');
        end
    else
        error('readGeodataBil:prismDate',['The PRISM date length of ' num2str(numel(numTest)) ' is not expected.']);
    end

    %Read test header:
    [hdrESRIRaw, metaHdr, bits] = read_PRISM_bin_hdr(fullfile(pathData,filePrismHdr{1}));
    for zz = 1 : numel(filePrismData)
        if all(~isnan([lonLd(:)', latLd(:)']))
            [dataCurr, hdrESRI, ~, ~] = read_PRISM_bin_data(fullfile(pathData,filePrismData{zz}), hdrESRIRaw, bits, [lonLd(:)', latLd(:)'], varReq, 'in');
        else
            [dataCurr, hdrESRI, ~, ~] = read_PRISM_bin_data(fullfile(pathData,filePrismData{zz}), hdrESRIRaw, bits, [], varReq, []);
        end

        if zz == 1
            sData.(varReq) = nan([numel(filePrismData), hdrESRI(2), hdrESRI(1)], 'single');
            dateLoad = nan([numel(filePrismData), 3]);
        end

        %Get date:
        indUnd = regexpi(filePrismData{zz}, '_');
        if regexpbl(prismType, 'clim')
            dateLoad(zz,2) = str2double(filePrismData{zz}(indUnd(end-1)+1 : indUnd(end)-1));
        else
            dateLoad(zz,1) = str2double(filePrismData{zz}(indUnd(end-1)+1 : indUnd(end-1)+4));
            dateLoad(zz,2) = str2double(filePrismData{zz}(indUnd(end-1)+5 : indUnd(end)-1));
        end

        sData.(varReq)(zz,:,:) = dataCurr;
    end


%         %If temperature, convert PRISM from Fahrenheit to Celsius:
%         if regexpbl(PrismVar, {'tmean', 'tmin', 'tmax'}) || regexpbl(varLd, {'tmp','tas', 'tmn', 'tmx'})
%             sData.(varLd) = (5/9)*(sData.(varLd) - 32);
%         end

    if regexpbl(prismType, 'clim')
        if ~issorted(dateLoad(:,2), 'rows')
            [dateLoad(:,2), indSort] = sortrows(dateLoad(:,2));
            sData.(varReq) = sData.(varReq)(indSort, :, :);
        end

        sData.(varReq) = squeeze(sData.(varReq));

        sData.info = {'source','PRISM'};
        sData.attTime = [sData.attTime; {'type', 'climatology'}];
    else
        dateRef = [1850,1,1];

        daysTest = days_since(dateRef, [dateLoad(:,1:2), 15*ones(numel(dateLoad(:,1)),1)], 'gregorian');
        if ~issorted(daysTest, 'rows')
            [~, indSort] = sort(daysTest);

            dateLoad = dateLoad(indSort,:);
            sData.(varReq) = sData.(varReq)(indSort, :, :);
        end

        sData.info = {'source','PRISM'};
        sData.attTime = [sData.attTime; {'type', 'time-series'}];
    end
    
%%WORLDCLIM
elseif strcmpi(type, 'WorldClim')
    %Check if this WorldClim dataset includes hyphen or not:
    if ~isempty(dir(fullfile(pathData, strcat(varLd, '_', '*.hdr' ) ) ))
        sepWC = '_';
    else
        sepWC = '';
    end

    %Find file for current month:
    if regexpbl(varLd,'alt')
        fileWcHdr = dir(fullfile(pathData, strcat(varLd, sepWC, '.hdr' ) ) );
        fileWcData = dir(fullfile(pathData, strcat(varLd, sepWC, '.bil' ) ) );
    else
        fileWcHdr = dir(fullfile(pathData, strcat(varLd, sepWC, num2str(mnthsLd), '.hdr' ) ) );
        fileWcData = dir(fullfile(pathData, strcat(varLd, sepWC, num2str(mnthsLd), '.bil' ) ) );
    end
    

    if isempty(fileWcHdr) || isempty(fileWcData) 
        error('readGeodataBil:BinMissing','The binary file was not found');
    elseif numel(fileWcHdr) > 1 || numel(fileWcData) > 1 
        error('readGeodataBil:BinExtra','Too many binary files were found.')
    end
    
    fileWcHdr  = extract_field(fileWcHdr, 'name');
    fileWcData = extract_field(fileWcData, 'name');
    
    for ii = 1 : numel(fileWcData(:))
        [hdrESRI, metaHdr] = read_wc_bin_hdr(fullfile(pathData,fileWcHdr{ii}));

        if all(~isnan([lonLd(:)', latLd(:)']))
            [dataTemp, hdrESRI, ~, ~] = read_wc_bin_data_v2(fullfile(pathData,fileWcData{ii}), hdrESRI, [lonLd(:)', latLd(:)'], varReq, 'in');
        else
            [dataTemp, ~, ~, ~] = read_wc_bin_data_v2(fullfile(pathData,fileWcData{ii}), hdrESRI, [], varReq, []);
        end
        
        if ii == 1
            sData.(varReq) = nan([numel(fileWcData(:)), size(dataTemp)], 'single');
            dateLoad = nan(numel(fileWcData(:)), 3);
        end
        
        sData.(varReq)(ii,:,:) = dataTemp;
        dateLoad(ii,2) = mnthsLd(ii);
    end

    
    sData.info = {'source','WorldClim'};
    sData.attTime = [sData.attTime; {'type', 'climatology'}];
    
%APHRODITE    
elseif strcmpi(type, 'APHRODITE')
    %Find files:
    fileAphroData = dir(fullfile(pathData, strcat('APHRO_', '*')));
    fileAphroData = extract_field(fileAphroData, 'name');
    
    %Get years
    yrsData = nan(numel(fileAphroData), 1);
    for ii = 1 : numel(fileAphroData)
        yrsData(ii) = str2double(fileAphroData{ii}(end-3:end));
    end
    indUse = find(yrsData >= min(yrsLd) & max(yrsLd));
    
    if numel(indUse) < max(yrsLd) - min(yrsLd) + 1
        warning('readGeodataBil:missingYrs',['The requested year range is ' ...
            num2str(min(yrsLd)) ' thru '  num2str(max(yrsLd)) ...
            ' but the years found are ' num2str(min(yrsData(indUse))) ...
            ' thru ' num2str(max(yrsData(indUse))) '.'])
    end
    
    fileAphroData = fileAphroData(indUse);
%     yrsData = yrsData(indUse);
    
    dateLoad = nan(0,3,'single');
    for ii = 1 : numel(fileAphroData)
        pathAphro = fullfile(pathData, fileAphroData{ii});
        [dataCurr, dateCurr, hdrESRI, metaHdr, ~, ~] = read_APHRODITE_bin(pathAphro, [lonLd(:)', latLd(:)'], mnthsLd, 'in');
        
        if ii == 1
            sData.(varReq) = nan([0, hdrESRI(2), hdrESRI(1)], 'single');
        end
        
        sData.(varReq) = cat(1, sData.(varReq), dataCurr);
        dateLoad = [dateLoad; dateCurr];
    end
    
    sData.info = {'source','APHRODITE'};
    sData.attTime = [sData.attTime; {'type', 'time-series'}];
end

sData.(varDate) = dateLoad;
[sData.(varLat), sData.(varLon)] = ESRI_hdr2geo(hdrESRI, metaHdr);
