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

function write_geodata_v2(pathData,sData,prec,wrtTyp,varargin)


%Parse variable input arguments:
blDisp = 1;
mnthsWrt = (1:12);
var = '';
if numel(varargin(:)) > 0
    for ii = 1 : numel(varargin(:))
        if regexpbl(varargin{ii}, {'no','disp'}, 'and')
            blDisp = 0;
        elseif regexpbl(varargin{ii}, 'months')
            mnthsWrt = varargin{ii+1};
        elseif regexpbl(varargin{ii}, 'var')
            var = varargin{ii+1};
        end
    end
end

varLon = 'longitude';
varLat = 'latitude';
varDate = 'date';

if blDisp == 1
    disp(['Data are now being written to ' pathData]);
end

%Create output directory if it doesn't already exist:
[foldIn, fileIn, ext] = fileparts(pathData);
if ~exist(foldIn,'dir')
    mkdir(foldIn);
end

if  ~isempty(regexpi(wrtTyp,'csv'))
    %Check that input data are single point:
    if sum(size(sData.(var)) ~= 1) > 1 || numel(sData.(varLat)) > 1 || numel(sData.(varLon)) > 1
        error('writeGeodata:csvNonPoint',['The data being written to ' ...
            pathData ' is not point data. It therefore cannot be written '...
            'to CSV file. Choose different output format.']);
    end
    
    ext = '.csv';
    [root,fileIn,~] = fileparts(pathData);
    pathData = fullfile(root,[fileIn, ext]);
    
    %Figure out time format:
    if isfield(sData, 'time')
        if sum(size(sData.time(1,:)) == 2) == 1 %time data likely stored as [year, month]
            vecTCurr = sData.time;
            if vecTCurr(1) < vecTCurr(2)
                vecTCurr(1:2) = fliplr(vecTCurr(1:2));
            end
        elseif sum(size(sData.time(1,:)) == 1) == 2 %Likely the time data is recorded as 'days since ...'
            [dateRef, tUnits] = NC_time_units(sData.atttime);

            if ~isempty(regexpi(tUnits,'days since'))
                cal =  NC_cal(sData.atttime);

                vecTCurr = days_2_date_v2(sData.time, dateRef, cal);
            else
                disp(['The GCM' char(39) 's time units are' tUnits '.']);
                error('NC_time:refUnknown',['A case has not been '...
                    'written to deal with the current time '...
                    'reference frame.']);
            end
        else
            error('writeGeodata:unknownTime',...
                ['The dimensions of the time vector are ' ...
                num2str(length(sData.time(:,1))) ' rows by ' ...
                num2str(length(sData.time(1,:))) ' columns.'  char(10) ...
                'A case to interpret this has not been written.'])
        end
    elseif isfield(sData, varDate)
        vecTCurr = sData.(varDate);
    else
        error('writeGeodata:noTimeFldCsv', 'No time/date field has been found.');
    end
    
    wrt_csv_ts(pathData, sData, vecTCurr, var);
    
elseif ~isempty(regexpi(wrtTyp,'asc'))
    ext = '.asc';
    [root,fileIn,~] = fileparts(pathData);
    pathData = fullfile(root,[fileIn,'.asc']);
    
    if isempty(prec)
        if regexpbl(char(sData.var),'pre')
            prec = 0;
        else
            prec = 1;
        end
    end
    
    if length(size(sData.(var))) > 2
        nIter = length(sData.(var)(:,1,1));
    else
        nIter = 1;
    end
    
    for ii = 1 : nIter
        if nIter > 1
            %Find current time:
            if isfield(sData, varDate) && numel(sData.(varDate)(:,1)) == nIter
                vecTCurr = sData.(varDate)(ii,:);
            elseif isfield(sData, 'time') && isfield(sData, 'attTime') 
                if sum(size(sData.time(ii,:)) == 2) == 1 %time data likely stored as [year, month]
                    vecTCurr = sData.time(ii,:);
                    if vecTCurr(1) < vecTCurr(2)
                        vecTCurr(1:2) = fliplr(vecTCurr(1:2));
                    end
                elseif sum(size(sData.time(ii,:)) == 1) == 2 %Likely the time data is recorded as 'days since ...'
                    [dateRef, tUnits] = NC_time_units(sData.atttime);

                    if ~isempty(regexpi(tUnits,'days since'))
                        cal =  NC_cal(sData.atttime);

                        vecTCurr = days_2_date(sData.time(ii), dateRef, cal);
                    else
                        disp(['The GCM' char(39) 's time units are' tUnits '.']);
                        error('NC_time:refUnknown',['A case has not been '...
                            'written to deal with the current time '...
                            'reference frame.']);
                    end
                else
                    error('write_geodata:unknownTime',...
                        ['The dimensions of the time vector are ' ...
                        num2str(length(sData.time(:,1))) ' rows by ' ...
                        num2str(length(sData.time(1,:))) ' columns.'  char(10) ...
                        'A case to interpret this has not been written.'])
                end
            else
                error('writeGeodata:noTimeFldAsc', 'No time/date field has been found.');
            end
            
            indDate = regexpi(fileIn,'_');
            if ~isempty(indDate)
                fileCurr = [fileIn, '_', num2str(vecTCurr(1)), ...
                    '_' num2str(vecTCurr(2)), ext];
            end
            pathCurr = fullfile(foldIn, fileCurr);
        else
            pathCurr = pathData;
        end
        
        %Check that step sizes are all equal.  
        step = nanmean(abs(horzcat(reshape(diff(sData.(varLon)),1,[]), reshape(diff(sData.(varLat)),1,[]))));
        %If stepsizes not equal, resample data:

        if round(100*(step - abs(diff(sData.(varLat)(1:2))))/step) ~= 0 || round(100*(step - abs(diff(sData.(varLon)(1:2))))/step) ~= 0
            %Display warning alerting about resampling 
            if ii == 1
                warning('write_geodata:avgStep',['The data currently being '...
                    'written have a non-uniform step size and will '...
                    'therefore be resampled to a uniform grid.  This will cause inaccurate '...
                    'representation.  To write data using non-uniform '...
                    'step-size, write in NetCDF format.']);
                
                %ADJUST LONGITUDE
                %Create longitude vector and parameters for:
                nLonPts = numel(sData.(varLon));
                stepAlign = 0.001*step;
                lonEdgO = box_edg(sData.(varLon));

                nLonAlignPts = round(2*step/stepAlign);
                lonLap = nan(nLonAlignPts,1);
                if sData.(varLon)(1) < sData.(varLon)(end)
                    for jj = 1 : nLonAlignPts
                        lonR1 = sData.(varLon)(1) - step + jj*stepAlign;
                        lonLap(jj) = sum(abs(lonEdgO - box_edg(lonR1 : step : lonR1 + (nLonPts-1)*step)));  
                    end
                elseif sData.(varLon)(1) > sData.(varLon)(end)
                    for jj = 1 : nLonAlignPts
                        lonR1 = sData.(varLon)(1) + step - jj*stepAlign;
                        lonLap(jj) = sum(abs(lonEdgO - box_edg(lonR1 : - step : lonR1 - (nLonPts-1)*step)));  
                    end
                else
                    error('write_geodata:lonOne','Case with unequal step size and only one indice not written for.')
                end

                %Find minimum displacement between original and new longitude
                %grids:
                [valLon, indLonLap] = min(lonLap);
                if valLon > nLonPts*0.5*step
                    warning('write_geodata:lonResample',['The total difference between longitude centroids values is ' num2str(valLon) ', which is larger than expected.']);
                end

                %Define new longitude grid:
                if sData.(varLon)(1) < sData.(varLon)(end)
                    lonR1 = sData.(varLon)(1) + step - indLonLap*stepAlign;
                    lonR = lonR1 : step : lonR1 + (nLonPts-1)*step;
                elseif sData.(varLon)(1) > sData.(varLon)(end)
                    lonR1 = sData.(varLon)(1) - step + indLonLap*stepAlign;
                    lonR = lonR1 : -step : lonR1 - (nLonPts-1)*step;
                else
                    error('write_geodata:lonOne','Case with unequal step size and only one indice not written for.');
                end

                %ADJUST LATITUDE
                %Create latitude vector:
                nLatPts = numel(sData.(varLat));
                latEdgO = box_edg(sData.(varLat));
                nLatAlignPts = round(2*step/stepAlign);
                latLap = nan(nLatAlignPts,1);

                if sData.(varLat)(1) < sData.(varLat)(end)
                    for jj = 1 : nLatAlignPts
                        latR1 = sData.(varLat)(1) - step + jj*stepAlign;
                        latLap(jj) = sum(abs(latEdgO - box_edg(latR1 : step : latR1 + (nLatPts-1)*step)'));  
                    end
                elseif sData.(varLat)(1) > sData.(varLat)(end)
                    for jj = 1 : nLatAlignPts
                        latR1 = sData.(varLat)(1) + step - jj*stepAlign;
                        latLap(jj) = sum(abs(latEdgO - box_edg(latR1 : -step : latR1 - (nLatPts-1)*step)'));  
                    end
                else
                    error('write_geodata:lonOne','Case with unequal step size and only one indice not written for.')
                end

                %Find minimum displacement between original and new longitude
                %grids:
                [valLat, indLatLap] = min(latLap);
                if valLat > nLatPts*0.5*step
                    warning('write_geodata:lonResample',['The total difference between longitude centroids values is ' num2str(valLat) ', which is larger than expected.']);
                end

                if sData.(varLat)(1) < sData.(varLat)(end)
                    latR1 = sData.(varLat)(1) - step + indLatLap*stepAlign;
                    latR = (latR1 : step : latR1 + (nLatPts-1)*step)';
                elseif sData.(varLat)(1) > sData.(varLat)(end)
                    latR1 = sData.(varLat)(1) + step - indLatLap*stepAlign;
                    latR = (latR1 : -step : latR1 - (nLatPts-1)*step)';

                else
                    error('write_geodata:latOne','Case with unequal step size and only one indice not written for.')
                end
            end
            
            
            
            %RESAMPLE DATA TO NEW GRID:
            if length(size(sData.(var))) == 3 && numel(latR) > 1 && numel(lonR) > 1
                dataCurr = geodata_area_wgt(sData.(varLon),sData.(varLat),sData.(var)(ii,:,:),lonR,latR);
            elseif length(size(sData.(var))) == 2 && numel(latR) > 1 && numel(lonR) > 1
                dataCurr = geodata_area_wgt(sData.(varLon),sData.(varLat),sData.(var)(:,:),lonR,latR);
            else
                dataCurr = NaN;
            end
            
            hdr = ESRI_hdr(lonR, latR, 'corner');
        else %Step sizes equal
            if length(size(sData.(var))) == 3
                dataCurr = squeeze(sData.(var)(ii,:,:));
            elseif length(size(sData.(var))) == 2
                dataCurr = sData.(var);
            else
                error('write_geodata:noData',['The number of grids found ' ...
                    'to write is ' num2str(nIter) '.  A case has not been '...
                    'written for this.']);
            end
            
            hdr = ESRI_hdr(sData.(varLon), sData.(varLat), 'corner');
        end

        write_ESRI_v4(squeeze(dataCurr), hdr, pathCurr, prec);
    end
elseif regexpbl(wrtTyp,{'netCDF','nc'})
    if isempty(ext) %If extension empty, create filename
        %Remove existing numbers from data path, if any, and replace with dates 
        %corresponding to all output:

        extNC = '.nc';
        indDate = regexpi(fileIn,'_');
        if ~isempty(indDate)
            fileIn = fileIn(1:indDate(end-1)-1);
        end

        %Determine if data are climatology:
        [blClim, yrsClim,~] = NC_isclim(sData);
        if blClim == 1
            yrsWrt = yrsClim;
            fileIn = [fileIn, '_clim'];
        else
            yrsWrt = [min(sData.(varDate)), max(sData.(varDate))];
        end
        
        %Ensure month string is two digits:
        if min(mnthsWrt) < 10
            mnthStart = ['0' num2str(min(mnthsWrt))];
        else
            mnthStart = num2str(min(mnthsWrt));
        end
        if max(mnthsWrt) < 10
            mnthEnd = ['0' num2str(max(mnthsWrt))];
        else
            mnthEnd = num2str(max(mnthsWrt));
        end
        
        %Find time step and use to get date for creating name
        %(Different date naming depending on monthly or daily data)
        if isfield(sData, 'date')
            if all(ismember(sData.(varDate)(:,2), mnthsWrt) ~= 0) && (numel(sData.(varDate)(1,:)) == 2 || all(sData.(varDate)(:,3) == sData.(varDate)(1,3)))
                timeStep = 'monthly';
            elseif numel(sData.('date')(1,:)) == 3 && mode(abs(diff(sData.date(:,3)))) == 1
                timeStep = 'daily';
            else
                timeStep = 'unknown';
            end
            
            if numel(sData.('date')(1,:)) == 2
                dayStart = 0;
                dayEnd = 0;
            else
                dayStart = sData.('date')(1,3);
                dayEnd = max(sData.('date')(:,3));
            end
        elseif isfield(sData, 'time')
            avgDiff = mode(diff(sort(sData.time)));

    %         avgDiff = nanmean(diff(sort(sData.time)));
            if avgDiff < 1.5 && avgDiff > 0.5
                timeStep = 'daily';
            elseif avgDiff < 45 && avgDiff > 10
                timeStep = 'monthly';
            elseif avgDiff < 380 && avgDiff > 350
                timeStep = 'yearly';
            elseif isnan(avgDiff)
               timeStep = 'unknown';
            else
                error('methodDirect:unknownTimeStep',['The current time step is ' ...
                    num2str(avgDiff) ', which has not been programmed for.']);
            end
            
            [dateRef, tUnits] = NC_time_units(sData.atttime);

            if ~isempty(regexpi(tUnits,'days since'))
                cal =  NC_cal(sData.atttime);

                vecTime = days_2_date(sData.time, dateRef, cal);
                
                dayStart = vecTime(1,3);
                dayEnd = max(vecTime(:,3));
            else
                disp(['The GCM' char(39) 's time units are' tUnits '.']);
                error('NC_time:refUnknown',['A case has not been '...
                    'written to deal with the current time '...
                    'reference frame.']);
            end
                    
        else
           timeStep = 'unknown';
        end
    
        %Create file name:
        if regexpbl(timeStep, 'daily')
            if dayStart < 10
                dayStart = ['0' num2str(dayStart)];
            else
                dayStart = num2str(dayStart);
            end
            if dayEnd < 10
                dayEnd = ['0' num2str(dayEnd)];
            else
                dayEnd = num2str(dayEnd);
            end
            
            fileIn = [fileIn, '_' num2str(yrsWrt(1)), mnthStart dayStart '-', ...
                num2str(yrsWrt(2)), mnthEnd, dayEnd, extNC];
        elseif regexpbl(timeStep, 'monthly')
            fileIn = [fileIn, '_' num2str(yrsWrt(1)), mnthStart '-', ...
                num2str(yrsWrt(2)), mnthEnd, extNC];
        elseif regexpbl(timeStep, 'unknown')
            fileIn = [fileIn '_unknown-dates'];
        else
            error('writeGeodata:unknownDates','The data appear to have an unexpected time resolution.')
        end

        %Full path:
        pathCurr = fullfile(foldIn, fileIn);
    else
        %Find time step
        if isfield(sData, 'date')
            if all(ismember(sData.(varDate)(:,2), mnthsWrt) ~= 0) && (numel(sData.(varDate)(1,:)) == 2 || all(sData.(varDate)(:,3) == sData.(varDate)(1,3)))
                timeStep = 'monthly';
            elseif numel(sData.('date')(1,:)) == 3 && mode(abs(diff(sData.date(:,3)))) == 1
                timeStep = 'daily';
            else
                timeStep = 'unknown';
            end
        elseif isfield(sData, 'time')
            avgDiff = mode(diff(sort(sData.time)));

    %         avgDiff = nanmean(diff(sort(sData.time)));
            if avgDiff < 1.5 && avgDiff > 0.5
                timeStep = 'daily';
            elseif avgDiff < 45 && avgDiff > 10
                timeStep = 'monthly';
            elseif avgDiff < 380 && avgDiff > 350
                timeStep = 'yearly';
            elseif isnan(avgDiff)
               timeStep = 'unknown';
            else
                error('methodDirect:unknownTimeStep',['The current time step is ' ...
                    num2str(avgDiff) ', which has not been programmed for.']);
            end      
        else
           timeStep = 'unknown';
        end
        
        pathCurr = pathData;
    end
    
    %Write data to NetCDF file:
    
    %Try to find data units:
    if isfield(sData, ['att' var])
        units = find_att(sData.(['att' var]), 'units', 'no_warning');
    else
        units = '';
    end
    
    %Create 'dateBnds':
    if regexpbl(timeStep, 'daily')
        dateBndsVec = sData.(varDate)(:,1:3);
    elseif regexpbl(timeStep, 'monthly')
        dateBndsVec = nan(2*numel(sData.(varDate)(:,1)), 3);
        cntr = 0;
        for ll = 1 : numel(sData.(varDate)(:,1))
            dateBndsVec(cntr+1:cntr+2,:) = [sData.(varDate)(ll,1), sData.(varDate)(ll,2), 1; ...
                sData.(varDate)(ll,1), sData.(varDate)(ll,2), eomday(sData.(varDate)(ll,1), sData.(varDate)(ll,2))];
            cntr = cntr + 2;
        end
    elseif regexpbl(timeStep, 'unknown')
        dateBndsVec = sData.(varDate)(:,1:3);
    end
    
    if isempty(var)
        error('write_geodata:noOutputVar', ['No output variable is present. '...
            'This can be provided via the optional input aregument and is '...
            'needed for writing to NetCDF files.']);
    end
    
    if ~isempty(units)
        print_grid_NC_v2(pathCurr, squeeze(sData.(var)), var, sData.(varLon), sData.(varLat), sData.(varDate), dateBndsVec, units, prec);
    else
        print_grid_NC_v2(pathCurr, squeeze(sData.(var)), var, sData.(varLon), sData.(varLat), sData.(varDate), dateBndsVec, prec);
    end
%     write_NC(pathCurr, sData, prec);
    
    %%THIS SECTION CASUES MATLAB TO CRASH EVERY TIME (DUE TO MEMORY?)
%     %If current time includes ending time-series element, ensure NetCDF
%     %data are in chronological order:
%     if sMeta.currTime(2) == sMeta.mnths(end)
%         %Read data from NetCDf file:
%         vecTime = ncread(pathCurr, 'time');
%         if ~issorted(vecTime)
%             keyboard
%             %Open NetCDF file and retrieve variable attributes:
%             fidNC = netcdf.open(pathCurr,'WRITE');
%             vidData = netcdf.inqVarID(fidNC,char(sData.var));
%             vifTime = netcdf.inqVarID(fidNC,'time');
%             aryData = netcdf.getVar(fidNC, vidData);
% 
%             %Sort time and data:
%             [vecTime, ordSrt] = sort(vecTime);
%             aryData = aryData(ordSrt,:,:);
% 
%             %Write re-ordered time and data vectors:
%             netcdf.putVar(fidNC,vifTime, 0, numel(vecTime), vecTime);
%             netcdf.putVar(fidNC,vidData, [0, 0, 0], size(aryData), aryData);
%             netcdf.close(fidNC);
%         end
%     end
else
   error('write_data:noType',['The data output format chosen is ' ...
       wrtTyp ', which is an unkown type.']); 
end


if blDisp == 1
    disp(['Data have finished being written to ' pathData]);
end