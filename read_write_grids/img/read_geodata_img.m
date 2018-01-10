function sData = read_geodata_img(pathData)

varLon = 'longitude';
varLat = 'latitude';

if ~regexpbl(pathData,'USGS')
   error('read_geodata:unknownIMG',['Currently, the only IMG ' ...
       'files that can be read are those from the USGS. ' char(10) ...
       'The reason is that it is known that the USGS file names ' ...
       'provide the position and resolution information.']);
end

%Find all files in folder:
if isdir(pathData)
    foldData = pathData;
else
    [foldData, ~, ~] = fileparts(pathData);
end


filesTemp = dir([foldData, filesep, '*', ext]);
filesTemp = extract_field(filesTemp,'name');

%LOAD ALL FILES IN FOLDER
I = cell(numel(filesTemp),1);
iLon = cell(numel(filesTemp),1);
iLat = cell(numel(filesTemp),1);

tempLon = [];
tempLat = [];
for ii = 1 : numel(filesTemp)
    indUnd = regexpi(filesTemp{ii}, '_');
    strCrd = filesTemp{ii}(indUnd(end-1)+1:indUnd(end)-1);

    if numel(strCrd) < 6
        strCrd = filesTemp{ii}(indUnd(end)+1:end-4);
    end

    res = filesTemp{ii}(indUnd(2)+1:indUnd(3)-1);
    if strcmpi(res, '2') %Image resolution is 2 arc seconds
        delta = 2/3600;
    else
        error('read_geodata:unknownImgRes',['The resolution ' res ...
            ' has not been programmed for.']);
    end

    indQuad = regexp(strCrd,'\D*');
    indOff = regexp(strCrd,'\d*');

    %Define latitude grid
    if strcmpi(strCrd(indQuad(1)), 'n')
        edgN = str2double(strCrd(indOff(1):indQuad(2)-1));
        edgS = edgN - 1; %1 degree tiles;
        iLat{ii} = (edgN - 0.5*delta : -delta : edgS + 0.5*delta)';
    elseif strcmpi(strCrd(indQuad(1)), 's')
        edgN = -str2double(strCrd(indOff(1):indQuad(2)-1));
        edgS = edgN - 1; %1 degree tiles;
        iLat{ii} = (edgN - 0.5*delta : -delta : edgS + 0.5*delta)';
    else
        error('read_geodata:unknownQuad',['The quadrant nomenclature ' ...
            strCrd(indQuad(1)) ' has not been programmed for.']);
    end

    %Define longitude grid:
    if strcmpi(strCrd(indQuad(2)), 'w')
        edgW = -str2double(strCrd(indOff(2):end));
        edgE = edgW + 1; %1 degree tiles;
        iLon{ii} = (edgW + 0.5*delta : delta : edgE - 0.5*delta);
    elseif strcmpi(strCrd(indQuad(2)), 'e')
        edgW = str2double(strCrd(indOff(2):end));
        edgE = edgW + 1; %1 degree tiles;
        iLon{ii} = (edgW + 0.5*delta : delta : edgE - 0.5*delta);
    else
        error('read_geodata:unknownQuad',['The quadrant nomenclature ' ...
            strCrd(indQuad(2)) ' has not been programmed for.']);
    end

    fid = fopen(fullfile(foldData, filesTemp{ii})); 
    I{ii} = fread(fid, [numel(iLat{ii}), numel(iLon{ii})]); 
    fclose(fid); 

    tempLon = [tempLon, iLon{ii}];
    tempLat = [tempLat; iLat{ii}];
end

%Find all unique lat and lon points
sData.(varLon) = unique(tempLon);
sData.(varLat) = flipud(unique(tempLat));

%Crop lat and lon pts to 
if isfield(sMeta,'crd') 
    sData.(varLon)(sData.(varLon) < sMeta.crd(1)) = [];
    sData.(varLon)(sData.(varLon) > sMeta.crd(2)) = [];
    sData.(varLat)(sData.(varLat) < sMeta.crd(3)) = [];
    sData.(varLat)(sData.(varLat) > sMeta.crd(4)) = [];
end

valNoData = -9999;
sData.data = valNoData*ones([numel(sData.(varLat)), numel(sData.(varLon))], 'like', I{1});

for ii = 1 : numel(I)
    %Check the ordering of the loaded geoTif:
    blLon = 0; %0 = ordered, 1 = flip, 2 = other
    blLat = 0; %0 = ordered, 1 = flip, 2 = other
    if ~issorted(iLon{ii})
        if issorted(fliplr(iLon{ii}))
            blLon = 1;
        else
            blLon = 2;
        end
    end
    if ~issorted(flipud(iLat{ii}))
        if issorted(iLat{ii})
            blLat = 1;
        else
            blLat = 2;
        end
    end


    %Find indices:
    indLonCurr = find(ismember(sData.(varLon), iLon{ii}));
    indLatCurr = find(ismember(sData.(varLat), iLat{ii}));

    if numel(indLonCurr) == numel(iLon{ii}) && numel(indLatCurr) == numel(iLat{ii})
        if blLat == 0 && blLon == 0
            sData.data(indLatCurr, indLonCurr) = I{ii};
        elseif blLat == 1 && blLon == 0
            sData.data(indLatCurr, indLonCurr) = flipud(I{ii});
        elseif blLat == 0 && blLon == 1
            sData.data(indLatCurr, indLonCurr) = fliplr(I{ii});
        elseif blLat == 1 && blLon == 1
            sData.data(indLatCurr, indLonCurr) = rot90(I{ii},2);
        elseif blLon == 2 || blLat == 2
            error('read_geodata:geoTifOrder','This geoTif indice ordering has not been programmed for.');
        end
    elseif ~isempty(indLonCurr) && ~isempty(indLatCurr) 
        [indLonCurr, indLonDataCurr] = ismember(sData.(varLon), iLon{ii});
        [indLatCurr, indLatDataCurr] = ismember(sData.(varLat), iLat{ii});

        if all(indLonCurr == 1 | indLonCurr == 0)
            indLonCurr = find(indLonCurr == 1);
            indLatCurr = find(indLatCurr == 1);
        end
        if all(indLonDataCurr == 1 | indLonDataCurr == 0)
            indLonDataCurr = find(indLonDataCurr == 1);
            indLatDataCurr = find(indLatDataCurr == 1);
        end

        if blLat == 0 && blLon == 0
            sData.data(indLatCurr, indLonCurr) = I{ii}(indLatDataCurr, indLonDataCurr);
        elseif blLat == 1 && blLon == 0
            sData.data(indLatCurr, indLonCurr) = flipud(I{ii}(indLatDataCurr, indLonDataCurr));
        elseif blLat == 0 && blLon == 1
            sData.data(indLatCurr, indLonCurr) = fliplr(I{ii}(indLatDataCurr, indLonDataCurr));
        elseif blLat == 1 && blLon == 1
            sData.data(indLatCurr, indLonCurr) = rot90(I{ii}(indLatDataCurr, indLonDataCurr),2);
        elseif blLon == 2 || blLat == 2
            error('read_geodata:geoTifOrder','This geoTif indice ordering has not been programmed for.');
        end
    end
end

sData.data(sData.data == valNoData) = nan;

outTyp = 'none';