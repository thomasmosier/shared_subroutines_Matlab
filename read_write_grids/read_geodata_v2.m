% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function sData = read_geodata_v2(pathData, varLd, lonLd, latLd, yrsLd, fr, cropType, varargin)


%Set select names of variables for output structure array:
varLat = 'latitude';
varLon = 'longitude';
varDate = 'date';


%List of extensions that can be read:
extRead = {'nc','asc','tif','img','bil','hdr'};

%Check if path is file or directory
[~,fileTemp,extTemp] = fileparts(pathData);
if exist(pathData, 'dir') ~= 0 && isempty(extTemp)
    files = extract_field(dir(pathData), 'name');
    ext = blanks(0);
    for ii = 1 : numel(files)
       if exist(fullfile(pathData, files{ii}),'file') && ~exist(fullfile(pathData, files{ii}),'dir') && ~strcmpi(files{ii}(1), '.')
           [~,~,ext] = fileparts(files{ii});
           break
       end
    end
    if ~regexpbl(ext, extRead)
        error('readGeodata:unknownExt', ['The input path ' char(39) fullfile(pathData,files{1}) ...
            char(39) ' does not have an extension that can be read by this program.']);
    end
elseif exist(pathData, 'file') ~= 0
    if isempty(extTemp) && regexpbl(pathData, 'gpcc')
        ext = 'gpcc';
    elseif ~isnan(str2double(extTemp))
        ext = 'numeric';
    else
        ext = extTemp;
    end
else 
    error('readGeodata:unknownPath', ['The input path ' char(39) pathData ...
        char(39) ' is not a directory or file.']);
end


%Parse variable input arguments:
blDisp = 1;
blMult = 1;
mnthsLd = (1:12);
units = 'same';
if numel(varargin(:)) > 0
    for ii = 1 : numel(varargin(:))
        if regexpbl(varargin{ii}, 'onefile')
           blMult = 0; 
        elseif regexpbl(varargin{ii}, {'no','disp'}, 'and')
            blDisp = 0;
        elseif regexpbl(varargin{ii}, 'months')
            mnthsLd = varargin{ii+1};
        elseif regexpbl(varargin{ii}, 'units')
            units = varargin{ii+1};
        end
    end
end

%Display message that file is about to be read
if blDisp == 1
    if all(~isnan(yrsLd))
        disp(['Data for months ' num2str(mnthsLd(1)) '-' num2str(mnthsLd(end)) ' are now being read from ' char(39) ...
            char(pathData) char(39) '.']);
    else
        disp(['Data are now being read from ' char(39) ...
            char(pathData) char(39) '.']);
    end
end
    
    
 
%READ NetCDF files:
if regexpbl(ext, 'nc')
    if blMult == 1
        sData = read_geodata_nc(pathData, lonLd, latLd, yrsLd, fr, cropType);
    else
        sData = read_geodata_nc(pathData, lonLd, latLd, yrsLd, fr, cropType, 'onefile');
    end

    varFld = fieldnames(sData);
    
    if isempty(varFld)
        warning('readGeodata:noFldsNc','No fields were loaded for current NetCDF file.');
        return
    end
    
    if any(strcmpi(varFld, varLd))
        varUse = varLd;
    else
        varUse = '';
        for ii = 1 : numel(varFld)
            if ~regexpbl(varFld{ii}, 'att')
                if strcmpi(varLd, 'pre')
                    if strcmpi(varFld{ii}, 'pr') || regexpbl(varFld{ii}, {'pre'})
                        varUse = varFld{ii};
                    end
                    
                    if isempty(varUse) && strcmpi(varFld{ii}, 'p')
                        varUse = varFld{ii};
                    end
                elseif strcmpi(varLd, 'tmp')
                    if strcmpi(varFld{ii}, 'tas') || regexpbl(varFld{ii}, {'tav','tmean'})
                        varUse = varFld{ii};
                    end  
                elseif strcmpi(varLd, 'tmn')
                    if strcmpi(varFld{ii}, 'tasmin') || regexpbl(varFld{ii}, {'tmin'})
                        varUse = varFld{ii};
                    end 
                elseif strcmpi(varLd, 'tmx')
                    if strcmpi(varFld{ii}, 'tasmax') || regexpbl(varFld{ii}, {'tmax'})
                        varUse = varFld{ii};
                    end
                elseif strcmpi(varLd, 'z')
                    if strcmpi(varFld{ii}, 'gp') || strcmpi(varFld{ii}, 'z') || regexpbl(varFld{ii}, 'dem') || regexpbl(varFld{ii}, 'orog')
                        varUse = varFld{ii};
                    end
                end
            end
        end
        clear ii

        if isempty(varUse)
           indUse = strcmpi(varFld, varLd);
           
            if ~any(indUse)
            	error('readGeodata:unknownVar',['The variable ' varLd ' has not been programmed for and was not found in the input file.']);
            else
                varUse = varLd;
            end
        end
    end
    
    if isempty(varUse)
        error('readGeodata:noNcVarFound',['No variable was found in the '...
            'current NetCDF file that corresponds to ' varLd '.']);
    end
    
    if ~isequal(varUse, varLd) && isfield(sData, ['att' varUse])
        sData.(['att' varLd]) = sData.(['att' varUse]);
        sData = rmfield(sData, ['att' varUse]);
    end 
        
       
    %Crop for requested months:
    if isfield(sData, varDate)
        indKeep = find(ismember(sData.(varDate)(:,2), mnthsLd) ==1);

        sData.(varLd) = sData.(varUse)(indKeep,:,:);
        if ~isequal(varLd, varUse)
           sData = rmfield(sData, varUse); 
        end
        sData.('time') = sData.('time')(indKeep);
        sData.(varDate) = sData.(varDate)(indKeep,:);
    end
    
    %Assign input variable name to requested variable name
    if ~isequal(varUse, varLd)
        fldsAll = fieldnames(sData);
        for ii = numel(fldsAll) : -1 : 1
            if regexpbl(fldsAll{ii}, varUse)
                fldNew = strrep(fldsAll{ii}, varUse, varLd);
                sData.(fldNew) = sData.(fldsAll{ii});
                sData = rmfield(sData, fldsAll{ii});
            end
        end
    end
%READ Multiple types of Tif files:    
elseif regexpbl(ext, 'tif')
    if regexpbl(pathData, {'wc','worldclim'})
        sData = read_worldclim_tif(pathData, varLd, mnthsLd, lonLd, latLd, fr, cropType);
        
        sData.info = {'source','WorldClim'};
        if isfield(sData, 'attTime')
            sData.attTime = [sData.attTime; {'type', 'climatology'}];
        else
            sData.attTime = {'type', 'climatology'};
        end
    elseif regexpbl(pathData,'srtm')
        sData = read_SRTM_tif(pathData, varLd, lonLd, latLd, fr, cropType);
        sData.(varDate) = nan(1,3);
        sData.('time') = nan;
        sData.(['att' varLd]) = {'units', 'm'};
    else
        sData = read_geotif(pathData, varLd, lonLd, latLd, fr, cropType);
        sData.(varDate) = nan(1,3);
        sData.('time') = nan;
    end
   
%READ WorldClim, PRISM, or APHRODITE binary (bil) files:
elseif ((regexpbl(ext,{'bil','hdr'}) && regexpbl(pathData, {'PRISM','worldclim','wc'})) || (strcmpi(ext,'numeric') && regexpbl(fileTemp, 'APHRO'))) && ~regexpbl(ext,'asc') %WorldClim binary
    sData = read_geodata_bil(pathData, varLd, lonLd, latLd, yrsLd, mnthsLd);
    
%READ USGS IMG files
elseif regexpbl(ext,{'img'})  
    sData = read_geodata_img(pathData);
    
    sData.(varLd) = sData.(varLd);
    sData = rmfield(sData, 'data');
    sData.(varDate) = nan(1,3);
    sData.('time') = nan;
    
%READ Willmott and Matsuura (i.e. University of Delaware gridded
%time-series) data
elseif regexpbl(pathData,{'will','maat','udel'}) || strcmpi(ext, 'numeric')
    if ~isdir(pathData)
        warning('read_geodata:WMfiles','The path for W&M data is a file, which is being converted to a directory.')
        pathData = fileparts(pathData);
    end

    [sData.(varLd), sData.(varDate), hdrESRI, metaESRI] = read_UDel_data(pathData, yrsLd, mnthsLd);
    sData.info = {'source', 'Willmott and Matsuura'};

    [sData.(varLat), sData.(varLon)] = ESRI_hdr2geo(hdrESRI, metaESRI);
    
    %Set time variable:
    sData.('attTime') = {'units', 'days since 1900-01-01'; ...
        'calendar', 'gregorian'};
    sData.time = days_since([1900,1,1], [sData.(varDate)(:,1:2), 15*ones([numel(sData.(varDate)(:,1)), 1])], 'gregorian');
    
%READ GPCC data (not in NetCDF format)    
elseif regexpbl(pathData,'gpcc') && strcmpi(ext, 'gpcc')
    if ~isdir(pathData)
%         warning('read_geodata:GPCCfiles','The path for the GPCC data is a file, which is being converted to a directory.')
        pathData = fileparts(pathData);
    end

        warning('off', 'GPCCFileFind:nWrong')
    [filesUse, sData.(varDate)] = GPCC_file_find(pathData, yrsLd, mnthsLd);
        warning('on', 'GPCCFileFind:nWrong')

    if numel(filesUse(:)) == 0
        error('readGeodata:noGpccFiles', ['No GPCC files were found at ' pathData]);
    end
    
    for ii = 1 : length(filesUse(:))
        fileCurr = fullfile(pathData,filesUse{ii});

        [dataCurr, hdrESRI, metaESRI, ~, ~] = read_GPCC_ASCII(fileCurr);

        if ii == 1
            sData.(varLd) = nan(length(filesUse(:)),length(dataCurr(:,1)),length(dataCurr(1,:)),'single'); 
        end
        sData.(varLd)(ii,:,:) = dataCurr;

    end
    sData.info = {'source','GPCC'};

    [sData.(varLat), sData.(varLon)] = ESRI_hdr2geo(hdrESRI, metaESRI);
    
    %Set time variable:
    sData.('attTime') = {'units', 'days since 1900-01-01'; ...
        'calendar', 'gregorian'};
    sData.time = days_since([1900,1,1], [sData.(varDate)(:,1:2), 15*ones([numel(sData.(varDate)(:,1)), 1])], 'gregorian');
    
    
%READ ESRI ASCII gridded files (from multiple sources, including PRISM)
elseif regexpbl(ext, 'asc')
    [dirData,~,ext]= fileparts(pathData);

    if isempty(ext)
        dirData = pathData;
    end
    
    if ~isempty(ext) && blMult == 0 %Path is file; load only this one
        [sData.(varLd), hdrESRI, metaESRI] = read_ESRI(pathData);
        
        %Check if data loaded has time embedded in filename
        %(*_year_month.asc)
        indNum = regexpi(pathData,'\d');
        if ~isempty(indNum)
            indTime = indNum(([0,diff(indNum)] ~= 1));
            indUnd = regexpi(pathData,'[_.]');
            if numel(indTime) > 1 && numel(indUnd) > 1
                %ASSUMES MONTHLY DATA BECAUSE ONLY LOOKS FOR TWO NUMBERS:
                sData.(varDate) = [str2double(pathData(indTime(end-1):indUnd(end - 1)-1)), str2double(pathData(indTime(end):indUnd(end)-1))];
            end
        end
    else %Find all times in directory provided
        if regexpbl(pathData, 'PRISM')
            if all(~isnan(yrsLd))
                [filesUse, sData.(varDate), hdrESRI, metaESRI] = PRISM_ASC_find(dirData, yrsLd, mnthsLd);
            else
                [filesUse, sData.(varDate), hdrESRI, metaESRI] = PRISM_ASC_find(dirData);
            end
        else
            %Find files:
            if all(~isnan(yrsLd))
                [filesUse, sData.(varDate), hdrESRI, metaESRI] = ASC_file_find(dirData, yrsLd, mnthsLd);
            else
                [filesUse, sData.(varDate), hdrESRI, metaESRI] = ASC_file_find(dirData);
            end
        end

       %Load data:
       if ~isempty(filesUse)
           for ii = 1 : length(filesUse)
                [dataCurr, ~, ~] = read_ESRI(filesUse{ii});

                if ii == 1 %Initialize data array
                    sData.(varLd) = nan(length(filesUse(:)),length(dataCurr(:,1)),length(dataCurr(1,:)));  
                end
                sData.(varLd)(ii,:,:) = dataCurr;
           end
       else
           sData.(varLd) = single.empty([0,0,0]);
           sData.(varLat) = single.empty(0);
           sData.(varLon) = single.empty(0);
           sData.time = single.empty(0);
           return
       end
    end

    if isfield(sData, 'date') && ~isfield(sData, 'time')
        %Set time variable:
        sData.('attTime') = {'units', 'days since 1900-01-01'; ...
            'calendar', 'gregorian'};
        sData.time = days_since([1900,1,1], [sData.(varDate)(:,1:2), 15*ones([numel(sData.(varDate)(:,1)), 1])], 'gregorian');
    end

    if regexpbl(pathData,'cru')
        sData.info = {'source','CRU'};
        warning('read_geodata:CRUScale',['The current data are in a folder with CRU.  ' ...
            char(10) 'CRU Ascii files are scaled by a factor of ten.  The current data are NOT being scaled.'])
%         sData.(varLd) = sData.(varLd)/10; %All CRU data scaled by factor of ten
    elseif regexpbl(pathData,'worldclim')
        sData.info = {'source','WorldClim'};
        if regexpbl(varLd,{'tmp','tas','tmin','tmax','tmn','tmx'})
            sData.(varLd) = sData.(varLd)/10; %WorldClim temperature data scaled by factor of ten
        end
    elseif regexpbl(pathData,'PRISM')
        sData.info = {'source','PRISM'};
    elseif regexpbl(pathData,'delta')
        sData.info = {'source','Downscaled data'};
    else
        sData.info = {'source','Unknown'};
    end
    
   [sData.(varLat), sData.(varLon)] = ESRI_hdr2geo(hdrESRI, metaESRI);
   
   if all(~isnan([lonLd(:); latLd(:)]))
        sMeta.crd = [min(lonLd), max(lonLd), max(latLd), min(latLd)];
        [~, indLatLd] = NC_lat_use_v4(sData.(varLat), {'units','degrees_North'}, sMeta);
        [~, indLonLd] = NC_lon_use_v4(sData.(varLon), {'units','degrees_East'}, sMeta);
       
        sData.(varLat) = sData.(varLat)(min(indLatLd):max(indLatLd));
        sData.(varLon) = sData.(varLon)(min(indLonLd):max(indLonLd));
        if ndims(sData.(varLd)) == 3
            sData.(varLd) = sData.(varLd)( : , min(indLatLd):max(indLatLd), min(indLonLd):max(indLonLd));
        elseif ndims(sData.(varLd)) == 2
            sData.(varLd) = sData.(varLd)(min(indLatLd):max(indLatLd), min(indLonLd):max(indLonLd));
        else
            error('readGeodata:ascdims',['The number of dimensions in the loaded ascii file is ' num2str(ndims(sData.(varLd))) ', but either 2 or 3 are expected.']);
        end
   end
   
else
    error('readGeodata:unknownType',['The following file is of an unknown type or format: ' pathData '.']);
end



%
if regexpbl(units, 'standard')
    if isfield(sData,['att' varLd])
        %Find timestep:
        if isfield(sData, 'date')
            if isempty(sData.date)
                error('readGeodata:emptyDate','The date vector is empty. This is likely because the file did not contain the requested dates.')
            end
            
            if all(ismember(sData.('date')(:,2), mnthsLd) ~= 0) && (numel(sData.('date')(1,:)) == 2 || sum(mode(sData.('date')(:,3)) == sData.('date')(:,3)) > 0.5*numel(sData.('time')))
                sData.timestep = 'monthly';
            elseif mode(abs(diff(sData.date(:,3)))) == 1 && (numel(sData.('date')(1,:)) == 3 || (numel(sData.('date')(1,:)) > 3 && mode(abs(diff(sData.date(:,4)))) == 0) )
                sData.timestep = 'daily';
            else
                sData.timestep = 'unknown';
            end
        elseif isfield(sData, 'time')
            if isempty(sData.time)
                error('readGeodata:emptyTime','The time vector is empty. This is likely because the file did not contain the requested dates.')
            end
            
            avgDiff = mode(diff(sort(sData.time)));

    %         avgDiff = nanmean(diff(sort(sData.time)));
            if avgDiff < 1.5 && avgDiff > 0.5
                sData.timestep = 'daily';
            elseif avgDiff < 45 && avgDiff > 10
                sData.timestep = 'monthly';
            elseif avgDiff < 380 && avgDiff > 350
                sData.timestep = 'yearly';
            elseif isnan(avgDiff)
               sData.timestep = 'unknown';
            else
                error('methodDirect:unknownTimeStep',['The current time step is ' ...
                    num2str(avgDiff) ', which has not been programmed for.']);
            end
        else
           sData.timestep = 'unknown';
        end
        
        sData = struct_2_standard_units(sData, varLd, sData.timestep);
%         if ~regexpbl(sData.timestep, 'unknown')
%             %Convert units:
%             sData = struct_2_standard_units(sData, varLd, sData.timestep);
%         else
%             warning('readGeodata:unknownUnitsUnknownTStep','Units are not being converted because the data time step is unknown');
%         end
    else
        warning('readGeodata:unknownUnitsNoAtt','Units are not being converted because no fields in data structure for attributes.');
    end
end

%Display message that file has been read
if blDisp == 1
    if all(~isnan(yrsLd))
        disp(['Data for months ' num2str(mnthsLd(1)) '-' num2str(mnthsLd(end)) ' have finished loading from ' char(39) ...
            char(pathData) char(39) '.']);
    else
        disp(['Data have finished loading from ' char(39) ...
            char(pathData) char(39) '.']);
    end
end

