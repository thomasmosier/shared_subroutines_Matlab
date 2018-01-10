function sNc = read_geodata_nc(pathNc, lonBnds, latBnds, yrBnds, fr, cropType, varargin)

%Make sure file exists:
if ~exist(pathNc, 'file')
    error('readGeodataNc:noPath',['The NetCDF file ' char(39) pathNc ...
        char(39) ' does not exist and therefore cannot be loaded.']);
end

%Variable input arguments
blMult = 1;
if numel(varargin(:)) > 0
    for ii = 1 : numel(varargin(:))
        if regexpbl(varargin{ii}, 'onefile')
           blMult = 0; 
        end
    end
end

%Turn off calendar warning if input data are APHRODITE
blAph = 0;
if regexpbl(pathNc, 'APHRO')
    blAph = 1;
    warning('off', 'ncCal:unknownCal');
end

[foldNc, ~, extNc] = fileparts(pathNc);

if isempty(extNc) && exist(pathNc, 'dir')
    disp('The input path is a directory. It is being searched for NetCDF files.');
    foldNc = pathNc;
    pathNc = extract_field(dir(fullfile(foldNc, '*.nc')), 'name');
    pathNc = fullfile(foldNc, pathNc{1});
    extNc = '.nc';
end

if ~regexpbl(extNc, 'nc')
    error('read_clim_mod_nc:unknownExt', ['The requested file' char(39) ...
        's extension is ' extNc ', but .nc is required.']);
end

filesAll = extract_field(dir(fullfile(foldNc, ['*' extNc])), 'name');

%Keep only files that have the same root as the given file:
indUnd = regexpi(filesAll{1}, '_');
if ~isempty(indUnd)
    fileNmRt = filesAll{1}(1:indUnd(end));
else
    fileNmRt = filesAll{1};
end

for ii = numel(filesAll(:)) : -1 : 1
    if ~regexpbl(filesAll{ii}, fileNmRt)
        filesAll(ii) = [];
    end
end

%Test type of NetCDf file (try to find function that can read time from name:
[yrsMod1, mnthsMod1, ~] = filename_time(pathNc);
    
if isempty(yrBnds)
    yrBnds = nan(1,2);
end

blLoaded = 0;
if numel(filesAll(:)) > 1 && blMult == 1
    if all2d(isnan(yrsMod1)) && all2d(isnan(mnthsMod1)) 
        warning('read_clim_mod_nc:unknownNcTime',['The time of the current '...
            'NetCDF file is unknown. This is due to the format of the file string.' char(10) ...
            'The function will not attempt to merge multiple files.' char(10) ...
            'Current file = ' pathNc]);
    elseif all(~isnan(yrBnds)) && (min2d(yrBnds) < yrsMod1(1) || max2d(yrBnds) > yrsMod1(2)) %Search for multiple files to merge
        %Choose which files to load and in which order (ensure data are
        %continuous)
        modUse = zeros(length(filesAll(:)),4);

        for ii = 1 : length(filesAll(:))
            modTimeTemp = filename_time(fullfile(foldNc, filesAll{ii}));
            modUse(ii,:) = [modTimeTemp(1,1), modTimeTemp(2,1), modTimeTemp(1,2), modTimeTemp(2,2)];

            %Do not load if ...
            if modUse(ii,2) < min2d(yrBnds) || modUse(ii,1) > max2d(yrBnds) %GCM ends before requested years or begins after requested years.
                modUse(ii,:) = NaN;
%             elseif modUse(ii,2) == min2d(yrBnds) %GCM ends in first year to load and does not contain requested months
%                 modUse(ii,:) = NaN;
%             elseif modUse(ii,1) == max2d(yrBnds) %GCM begins in last year to load and does not contain requested months
%                 modUse(ii,:) = NaN;
            end
        end

        %Delete files that do not contain necessary time-series elements:
        filesAll(isnan(modUse(:,1))) = [];
        modUse(isnan(modUse(:,1)),:) = [];

        %Find temporal ordering:
        [~, ordStrtGcm] = sort(modUse(:,1));
        [~, ordEndGcm] = sort(modUse(:,2));
        if ~isequal(ordStrtGcm,ordEndGcm)
            warning('read_geodata:GcmOrder',['The temporal ordering of '...
                'GCMs to load depends on whether the start date or end '...
                'dates are compared.  The ordering using starts dates '...
                'will be used.'])
        end

        %Order files based upon dates:
        filesAll = filesAll(ordStrtGcm);
        modUse = modUse(ordStrtGcm,:);
        
        if any(modUse(2:end,1) ~= modUse(1:end-1,2)+1)
           warning('read_geodata:dataGap', ['There appears to be a gap '...
               'in the NetCDF data being loaded from ' foldNc '.']); 
        end

        sNc = load_and_merge_nc(fullfile(foldNc, filesAll), lonBnds, latBnds, yrBnds, fr, cropType);
        blLoaded = 1;
    end
end

if blLoaded == 0
    sNc = read_grid_NC_w_bnds_v2(pathNc, lonBnds, latBnds, yrBnds, fr, cropType);
end


%Make sure time/dates are sorted in chronological order:
if isfield(sNc, 'time') && ~issorted(sNc.time)
   [~, srtTime] = sort(sNc.time);

   fldsChk = fields(sNc);
   nTime = numel(sNc.time);
   for ii = 1 : numel(fldsChk(:))
       szCurr = size(sNc.(fldsChk{ii}));
       indTime = find(szCurr == nTime);
       if ~isempty(indTime)
           if numel(indTime) == 1
               if numel(szCurr) == 1
                   switch indTime
                        case 1
                            sNc.(fldsChk{ii}) = sNc.(fldsChk{ii})(srtTime);
                       otherwise
                           error('readGeodataNc:sz1IndTimeLarge',['The time dimension ' ...
                               'should be 1 or less, but it is ' num2str(indTime) '.']);
                   end
               elseif numel(szCurr) == 2
                    switch indTime
                        case 1
                            sNc.(fldsChk{ii}) = sNc.(fldsChk{ii})(srtTime,:);
                        case 2
                            sNc.(fldsChk{ii}) = sNc.(fldsChk{ii})(:,srtTime);
                        otherwise
                            error('readGeodataNc:sz2IndTimeLarge',['The time dimension '...
                                'should be 2 or less, but it is ' num2str(indTime) '.']);
                    end
               elseif numel(szCurr) == 3
                    switch indTime
                        case 1
                            sNc.(fldsChk{ii}) = sNc.(fldsChk{ii})(srtTime,:,:);
                        case 2
                            sNc.(fldsChk{ii}) = sNc.(fldsChk{ii})(:,srtTime,:);
                        case 3
                            sNc.(fldsChk{ii}) = sNc.(fldsChk{ii})(:,:,srtTime);
                        otherwise
                            error('readGeodataNc:sz3IndTimeLarge',['The time dimension '...
                                'should be 3 or less, but it is ' num2str(indTime) '.']);
                    end
               else
                   error('readGeodataNc:szCurrLarge',['The current field, ' fldsChk{ii} ...
                       ', has ' num2str(numel(szCurr)) ' dimensions, which has not been programmed for.']);
               end
           else
              error('readGeodataNc:indTimeUnknown',['More than one dimension is field ' ...
                  fldsChk{ii} ' has the same number of elements as the time vector. ' ...
                  'This case is unexpected and has not been programmed for.']);
           end
       end
   end
   clear ii
end

%Change 'lon' to 'longitude' and 'lat' to 'latitude'
if isfield(sNc, 'lon')
   sNc.longitude = sNc.lon;
   sNc = rmfield(sNc, 'lon');
end
if isfield(sNc, 'Lon')
   sNc.longitude = sNc.Lon;
   sNc = rmfield(sNc, 'Lon');
end
if isfield(sNc, 'Longitude')
   sNc.longitude = sNc.Longitude;
   sNc = rmfield(sNc, 'Longitude');
end
if isfield(sNc, 'lat')
   sNc.latitude = sNc.lat;
   sNc = rmfield(sNc, 'lat');
end
if isfield(sNc, 'Lat')
   sNc.latitude = sNc.Lat;
   sNc = rmfield(sNc, 'Lat');
end
if isfield(sNc, 'Latitude')
   sNc.latitude = sNc.Latitude;
   sNc = rmfield(sNc, 'Latitude');
end


%Turn warning back on (if turned off previously
if blAph == 1
    warning('on', 'ncCal:unknownCal');
end

% %Check that requested data field present
% if ~isfield(sEra, varData)
%    error('read_ERAi_nc:missingVar',['The requested data variable, ' ...
%        varData ', is not present in the loaded file.']); 
% end

% %Reorder dimensions (put time first)
% if ndims(sEra.(varData)) == 3
%     sEra.(varData) = permute(sEra.(varData), [3,2,1]); 
%     indFlip = 2;
% elseif ismatrix(sEra.(varData))
%     sEra.(varData) = permute(sEra.(varData), [2,1]); 
%     indFlip = 1;
% else
%     warning('read_ERAi_nc:unknownDim', ...
%         ['The input array dimensions are ' ...
%         num2str(ndims(sEra.(varData))) ', which is unexpected.'])
% end

%Previously it appeared that latitude of ERA arrays was inverted...
% sEra.(varData) = flip(sEra.(varData), indFlip); 
% 
% %Output longitude as column vector
% sEra.longitude = sEra.longitude(:)';

% 
% %Calculate positions in UTM projection:
% [sEra.(varUtmN), sEra.(varUtmE), sEra.(varUtm)] = deg2utm_grid(sEra.latitude, sEra.longitude);

%figure; pcolor(sEra.(varData)); colorbar
% figure; pcolor(sEra.longitude, sEra.latitude, sEra.(varData)); colorbar