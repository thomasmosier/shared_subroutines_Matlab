function sData = read_grid_NC_w_bnds_v2(path, lonBnds, latBnds, yrBnds, fr, cropType)


% varRem = {'lon', 'longitude', 'lat', 'latitude', 'time', 'time_bnds'};
% vatLat = {'lat', 'latitude'};
% vatLon = {'lon', 'longitude'};

% %Find possible data variables:
% if isempty(varargin) && isempty(varargin{1})
%     fInfo = ncinfo(path);
%     varLd = extractfield(fInfo, 'Variables');
%     if iscell(varLd) && isstruct(varLd{1})
%        varLd = extractfield(varLd{1},'Name'); 
%     end
% 
%     for ii = 1 : numel(varRem)
%         varLd(strcmpi(varLd, varRem{ii})) = [];
%     end
% else
%     vadLd = varargin{1};
% end

%Initialize output:
sData = struct;

%Check that path exists:
if ~exist(path,'file')
    error('read_grid_NC_w_bnds:pathNoExist',['The requested path does not exist: ' path]);
end

%Load general attributes:
attGen = ncinfo(path);
sData.('attributes') = squeeze(struct2cell(attGen.Attributes))';


%Load variables:
fInfo = ncinfo(path);
varLd = extract_field(fInfo.Variables, 'Name');
% varLd = extractfield(fInfo.Variables, 'Name');


%Find latitude variable:
varLatTest = {'lat','latitude','row','y','Lat','Latitude','Row','Y','LAT','LATITUDE','ROW'};
varLat = blanks(0);
for ii = 1 : numel(varLd)
    if any(strcmpi(varLatTest, varLd{ii}))
        varLat = varLd{ii};
        varLd(ii) = []; 
        break
    end
end

if isempty(varLat)
   error('readGridNcWBnds:noLat',['No latitude variable was found in the NetCDF file being loaded: ' path]) 
end

%Find longitude variable:
varLonTest = {'lon','longitude','col','x','Lon','Longitude','Col','X','LON','LONGITUDE','COL'};
varLon = blanks(0);
for ii = 1 : numel(varLd)
    if any(strcmpi(varLonTest, varLd{ii}))
        varLon = varLd{ii};
        varLd(ii) = []; 
        break
    end
end

if isempty(varLon)
   error('readGridNcWBnds:noLon',['No longitude variable was found in the NetCDF file being loaded: ' path]) 
end

%Load longitude, latitude, and time coordinates:
lonIn = ncread(path, varLon);
latIn = ncread(path, varLat);

%Ensure consistent orientation of lat and lon:
if any(size(lonIn) == 1) && any(size(latIn) == 1)
    lonIn = lonIn(:)';
    latIn = latIn(:);
    
    %Get sizes:
    nLonOrg = numel(lonIn);
    nLatOrg = numel(latIn);

    if ~all(diff(lonIn) > 0) && ~all(diff(lonIn) < 0) 
        error('readGridNcWBnds:outOfOrderLon','The longitude indices are out of order. This requires further programming.');
    end
    if ~all(diff(latIn) > 0) && ~all(diff(latIn) < 0) 
        error('readGridNcWBnds:outOfOrderLat','The latitude indices are out of order. This requires further programming.');
    end
    if all(size(lonIn) > 1) 
        error('readGridNcWBnds:szLon','The longitude field is a 2D array. This requires further programming.');
    end
    if all(size(latIn) > 1) 
        error('readGridNcWBnds:szLat','The latitude field is a 2D array. This requires further programming.');
    end
else
    %Ensure lat and lon are ordered as expected:
    lonRaw = extract_field(ncinfo(path, varLon), 'Dimensions');
    lonDim = extract_field(lonRaw{1}, 'Name');
    latRaw = extract_field(ncinfo(path, varLat), 'Dimensions');
    latDim = extract_field(latRaw{1}, 'Name');

    varLonTest = {'lon','longitude','col','x','Lon','Longitude','Col','X','LON','LONGITUDE','COL','west_east','east_west'};
    varLatTest = {'lat','latitude','row','y','Lat','Latitude','Row','Y','LAT','LATITUDE','ROW','north_south','south_north'};
        
    if ~regexpbl(lonDim{1}, varLatTest) && regexpbl(lonDim{1}, varLonTest)
        lonIn = permute(lonIn, [2,1]);
    elseif ~regexpbl(lonDim{1}, varLatTest) && ~regexpbl(lonDim{1}, varLonTest)
        error('readGridNcWBnds:unknownLonOrder', [lonDim{1} ' was not recognized as a latitude dimension.']);
    end
        
    if ~regexpbl(latDim{1}, varLatTest) && regexpbl(latDim{1}, varLonTest)
        latIn = permute(latIn,[2,1]);
    elseif ~regexpbl(latDim{1}, varLatTest) && ~regexpbl(latDim{1}, varLonTest)
        error('readGridNcWBnds:unknownLatOrder', [latDim{1} ' was not recognized as a latitude dimension.']);
    end
    
    nLonOrg = numel(lonIn(1,:));
    nLatOrg = numel(latIn(:,1));
end



%Create 'sMeta' struct, used in lat/lon indice functions
sMeta.crd = [sort(lonBnds(:)'), sort(latBnds(:)')];
if all(isnan(sMeta.crd))
    sMeta.crd = nan(1,4);
end
sMeta.frame = fr;
sMeta.cropType = cropType;

%Find lat/lon indices to load
attLon = ncinfo(path, varLon);
sData.(['att' varLon]) = squeeze(struct2cell(attLon.Attributes))';
attLat = ncinfo(path, varLat);
sData.(['att' varLat]) = squeeze(struct2cell(attLat.Attributes))';

if any(size(lonIn) == 1) && any(size(latIn) == 1)
    [latTest, indLatLd] = NC_lat_use_v4(latIn, sData.(['att' varLat]), sMeta);
    [lonTest, indLonLd] = NC_lon_use_v4(lonIn, sData.(['att' varLon]), sMeta);

    if latIn(indLatLd(1)) < latIn(indLatLd(end)) && indLatLd(1) < indLatLd(end)
        indLatLd = fliplr(indLatLd(:)');
        latTest = fliplr(latTest(:)');
    end
    
    if lonIn(indLonLd(1)) > lonIn(indLonLd(end)) && indLonLd(1) < indLonLd(end)
        indLonLd = fliplr(indLonLd(:)');
        lonTest = fliplr(lonTest(:)');
    end
else
    warning('readGridNc:matrixCrd', ['The lat or lon variable for the current NetCDF file is a matrix instead of a vector.' ...
        char(10) 'Cropping is not available for this case.']);
    
    latBnds = nan(1,2);
    lonBnds = nan(1,2);

    indLatLd = (1:numel(latIn(:,1)));
    indLonLd = (1:numel(lonIn(1,:)));
    
    if (min(latIn(1,:)) < max(latIn(end,:))) && indLatLd(1) < indLatLd(end)
        indLatLd = fliplr(indLatLd(:)');
    end

    if (min(lonIn(:,1)) > max(lonIn(:,end))) && indLonLd(1) < indLonLd(end)
        indLonLd = fliplr(indLonLd(:)');
    end
end



%Make sure lat and lon are in map orientation. Rotate these fields and
%their corresponding indices:
% if min(lonIn(:,1)) > max(lonIn(:,end))
% % 	if all(isnan(lonBnds))
% %         indLonLd = fliplr(indLonLd);
% % 	end
%     lonIn = fliplr(lonIn);
%     latIn = fliplr(latIn);
% end
% if min(latIn(1,:)) < max(latIn(end,:))
% % 	if all(isnan(latBnds))
% %         indLatLd = fliplr(indLatLd);
% % 	end
%     lonIn = flipud(lonIn);
%     latIn = flipud(latIn);
% end


%Crop lon:
sData.(varLon) = lonIn(:,indLonLd);
%If longitude goes from 0 -> 360, transform to -180 -> 180
indGT180 = find(sData.(varLon) > 180);
if ~isempty(indGT180)
    sData.(varLon)(indGT180) = sData.(varLon)(indGT180) - 360;
end
indEq180 = find(sData.(varLon) >= 180);
if ~isempty(indEq180) && indEq180(1) == 1 && diff(sData.(varLon)(1:2)) < 0
    sData.(varLon)(1) = sData.(varLon)(1) - 360;
end
if ~isempty(indEq180) && numel(indEq180) > 2
    error('readGridNcWBnds:GT180','More than two lon indices are 180. This is unexpected in regular grids.');
end

if all(~isnan(lonBnds)) && ~all(lonTest(:)' == sData.(varLon))
    error('readGridNcWBnds:lonMismatch','There is disagreement between the test and actual lon grid.');
end

%Crop lat:
sData.(varLat) = latIn(indLatLd,:);
    
if all(~isnan(latBnds)) && ~all(latTest(:) == sData.(varLat))
    error('readGridNcWBnds:latMismatch','There is disagreement between the test and actual lat grid.');
end



% %Find indices to load based on lat/lon variables:
% [indLonLd, indLatLd] = find_crop_ind(sData.(varLon), sData.(varLat), lonBnds, latBnds, fr, cropType);


%Find time variable:
varT = blanks(0);
for ii = 1 : numel(varLd)
    if regexpbl(varLd{ii}, {'time'})
        varT = varLd{ii};
        varLd(ii) = []; 
        break
    end
end


%Create date vector for time:
%Load time and select time indices to keep
if ~isempty(varT)
    sData.(varT) = double(ncread(path, varT));

    [sData.(['att' varT]), ~] = load_ncatt(path, varT);

    %Read time units from GCM file:
    [gcmRef, tUnits] = NC_time_units(sData.(['att' varT]));

    if regexpbl(gcmRef,'unknown')
        sData.date = nan(1,3);
    else
        %Find calendar info:
        strCal =  NC_cal(sData.(['att' varT]));

        if regexpbl(strCal, 'unknown')
            strCal = 'gregorian';
            if ~regexpbl(path, 'aphro')
                warning('readGridNcWBnds:unknownCal','The calendar is unknown, but is being set to Gregorian.');
            end
        end
        
        %Find time offset from the reference date:
        if regexpbl(tUnits, 'days since')
            sData.date = days_2_date_v2(sData.(varT), gcmRef, strCal);
        elseif regexpbl(tUnits, 'Ka BP') && regexpbl(strCal,'noleap')
            sData.date = time_kaBP_2_standard(sData.(varT),gcmRef,strCal);
        elseif regexpbl(tUnits, 'hours since')
            if numel(gcmRef) == 3
                gcmRef = [gcmRef, 0];
            end

            sData.date = days_2_date_v2(sData.(varT)/24, gcmRef, strCal);
        else
            warning('NC_time:unitsUnknown', ['The NetCDF file' char(39) ...
                's time units are' tUnits '. This has not been programmed for. '...
                'Therefore, time cropping is being disabled.']);
            
            yrBnds = nan(1,2);
            if regexpbl(tUnits, 'Model time in UTC, format YYYY-MM-DD_HH:MM')
                sData.(varT) = sData.(varT)(1,:);
                sData.(varT) = sData.(varT)(:);
            end
            
            sData.date = nan([numel(sData.(varT)), 3]);
        end
    end

    %Find time indices for cropping inputs
    if all(~isnan(yrBnds))
        indTime = find(sData.date(:,1) >= min(yrBnds) & sData.date(:,1) <= max(yrBnds));
        if isempty(indTime)
            sData.(varT) = [];
            sData.date = nan(0, numel(sData.date(1,:)));
%             warning('read_clim_mod_nc:noTime','No time indices within the requested yers have been found.') 
        end
    else
        indTime = (1 : numel(sData.(varT)));
    end
    
    %Crop time indices
    sData.(varT) = sData.(varT)(indTime,:);
    sData.date = sData.date(indTime,:);
else
    indTime = nan(2,1);
    warning('readGridNcWBnds:noTime','No time variable was found in the NetCDF file being loaded.');
end


for ii = 1 : numel(varLd)
    infoCurr = ncinfo(path, varLd{ii});
    szCurr = extract_field(infoCurr, 'Size');
%     if iscell(szCurr) 
%         if numel(szCurr) == 1
%             szCurr = szCurr{1}; 
%         else
%             error('readGridNcBnds:sizeCell','The size of the current field is a cell with more than one entry.')
%         end
%     end
    if isempty(szCurr{1})
        szCurr = extract_field(infoCurr, 'size');
    end
    if isempty(szCurr{1})
        warning('readGridNcBnds:noVarSize','Size attribute not found. All data will be loaded.');
    else
        szCurr = szCurr{1};
        
        if numel(szCurr) == 3
            sData.(varLd{ii}) = NC_field_load(path, varLd{ii}, indTime, indLatLd, indLonLd);
        elseif numel(szCurr) == 2
            %Find order of dimensions:
            dLat = nan;
            dLon = nan;
            if any(szCurr == nLonOrg) && any(szCurr == nLatOrg) 
                for jj = 1 : numel(szCurr)
                    switch szCurr(jj)
                        case nLatOrg
                            dLat = jj;
                        case nLonOrg
                            dLon = jj;
                        otherwise
                            warning('readGridNcBnds:noSizeDim',['There is no '...
                                'matching size dimension for ' num2str(jj) '.']);
                    end
                end

                %Find ordering:
                if numel(szCurr) == 2 && all(~isnan([dLat, dLon]))
                    [~, ordr] = sort([dLat, dLon]);

                    start = [min(indLatLd), min(indLonLd)];
                    count = [max(indLatLd) - min(indLatLd) + 1, max(indLonLd) - min(indLonLd) + 1];

                    %Load segment:
                    sData.(varLd{ii}) = ncread(path, varLd{ii}, start(ordr), count(ordr));

                    %Permute loaded data to desired orientation:
                    sData.(varLd{ii}) = permute(sData.(varLd{ii}), ordr); 
                    
                    %Check sorting of lat and lon
                    if ~issorted(indLatLd)
                        [~, indLatSrt] = sort(indLatLd);
                        
                        sData.(varLd{ii}) = sData.(varLd{ii})(indLatSrt,:);
                    end
                    if ~issorted(indLonLd)
                        [~, indLonSrt] = sort(indLonLd);
                        
                        sData.(varLd{ii}) = sData.(varLd{ii})(:,indLonSrt);
                    end
                end
            elseif any(szCurr == nLonOrg) 
                dimLon = find(szCurr == nLonOrg);
                
                if dimLon == 1
                    start = [min(indLonLd), 1];
                    count = [max(indLonLd) - min(indLonLd) + 1, szCurr(2)];
                elseif dimLon == 2
                    start = [1, min(indLonLd)];
                    count = [szCurr(1), max(indLonLd) - min(indLonLd) + 1];
                end
                
                temp = ncread(path, varLd{ii}, start, count);
                
                if ~issorted(indLonLd)
                    [~, indLonSrt] = sort(indLonLd);
                    if dimLon == 1
                        sData.(varLd{ii}) = temp(indLonSrt,:);
                    elseif dimLon == 2
                        sData.(varLd{ii}) = temp(:,indLonSrt);
                    end
                else
                    sData.(varLd{ii}) = temp;
                end
            elseif any(szCurr == nLatOrg)
                dimLat = find(szCurr == nLatOrg);
                
                if dimLat == 1
                    start = [min(indLatLd), 1];
                    count = [max(indLatLd) - min(indLatLd) + 1, szCurr(2)];
                elseif dimLat == 2
                    start = [1, min(indLatLd)];
                    count = [szCurr(1), max(indLatLd) - min(indLatLd) + 1];
                end
                
                temp = ncread(path, varLd{ii}, start, count);
                
                if ~issorted(indLatLd)
                    [~, indLatSrt] = sort(indLatLd);
                    
                    if dimLat == 1
                        sData.(varLd{ii}) = temp(indLatSrt,:);
                    elseif dimLat == 2
                        sData.(varLd{ii}) = temp(:,indLatSrt);
                    end
                else
                    sData.(varLd{ii}) = temp;
                end
            else
                sData.(varLd{ii}) = ncread(path, varLd{ii});
            end
        else
            sData.(varLd{ii}) = ncread(path, varLd{ii});
        end
    end
    
    sData = add_ncatt_2_struct(sData, path, varLd{ii});
end




% for jj = 1 : varLon
%     if regexpbl(varLon{jj}, varLd)
%         
%     end
% end
% lonMODIS = ncread(path, 'lon');
%     lonMODIS = lonMODIS(:)';
% latMODIS = ncread(path, 'lat');
%     latMODIS = latMODIS(:);
% MODISTemp = ncread(path, varLd);

