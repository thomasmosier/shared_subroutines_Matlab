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

function pathSnow = snowfall(foldPre,foldTmp, snowMod, param, fileTyp,wrtTs,varargin)
%Script finds info from ASCII files in selected folders and calculates snowfall
%INPUTS:
%fileTyp: ('ascii' or 'netcdf') is input file format.
%wrtTs: [0 or 1] and specifies if snowfall time-series grids should be written.
%'varargin': If netcdf selected, varargin must provide yrs and bounding
%coordinates to select.


%Initialize output
pathSnow = [];

%Confirm type of snow model
if regexpbl(snowMod,'step')
    if numel(param) ~= 1
       error('snowfall:nStep','For the step function snowfall mode, the parameter input variable must have length = 1.'); 
    end
elseif regexpbl(snowMod,'ramp')
    if numel(param) ~= 2
       error('snowfall:nRamp','For the linear ramp snowfall mode, the parameter input variable must have length = 2.'); 
    end
    
    if ~issorted(param)
       param = param([2,1]);
       warning('snowfall:prmFlip','Snow model parameters being flipped so that T_{snow} < T_{rain}.');
    end
else
   error('snowfall:snowMod',[char(39) snowMod char(39) ' is not a known snow model type.']) 
end

%Confirm if monthly output grids to be written:
if ~regexpbl({'asc','cdf'},fileTyp)
	error('snowfall:wrtTyp',[char(39) fileTyp char(39) ' was selected, '...
      'which is not an acceptable output file type.']) 
end



%If varargin for yrs and bounding coordinates, retrieve values:
if ~isempty(varargin) && ~isempty(varargin{1})
    yrs = [];
    crd = [];
    for ii = 1 : length(varargin(:))
       if length(varargin{ii}) == 2 %Argument is yrs
           yrs = varargin{ii};
       elseif length(varargin{ii}) == 4 %Argument is coordinates
           crd =  varargin{ii};
       end
    end
    if isempty(yrs) || isempty(crd)
        error('snowfall:varargin', ['Varargin are needed to specify ' ...
            'yrs and coordinates.  At least one of those variables is empty.'])
    end
end

%Display message confirming start of processing:
disp(['Snowfall data processing has begun.' char(10) ...
    'Precipitation data will be read from ' char(39) foldPre char(39) '.' char(10) ...
    'Mean Temperature data will be read from ' char(39) foldTmp char(39) '.']);

%Find files in directories:
if ~isempty(regexpi(fileTyp,'asc'))
    extTyp = 'asc';
elseif ~isempty(regexpi(fileTyp,'cdf'))
    extTyp = 'nc';
end


%Find files in each folder:
gridTmpFiles = find_files(foldTmp,extTyp);
gridPreFiles = find_files(foldPre,extTyp);

%%LOAD FILES AND REQUIRED ASSOCIATED OTHER DATA:
if regexpbl(fileTyp,'asc')    
    %Find Dates of each file:
    dateTmp = nan(numel(gridTmpFiles),2);
    datePre = nan(numel(gridPreFiles),2);
    for ii = 1 : numel(gridTmpFiles)
        indUnd = regexpi(gridTmpFiles{ii},'_');
        indExt = regexpi(gridTmpFiles{ii},'\.');
        dateTmp(ii,1) = str2double(gridTmpFiles{ii}(indUnd(end-1)+1:indUnd(end)-1));
        if isnan(dateTmp(ii,1))
            dateTmp(ii,1) = -9999;
        end
        dateTmp(ii,2) = str2double(gridTmpFiles{ii}(indUnd(end  )+1:indExt(end)-1)); 
    end
    for ii = 1 : numel(gridPreFiles)
        indUnd = regexpi(gridPreFiles{ii},'_');
        indExt = regexpi(gridPreFiles{ii},'\.');
        datePre(ii,1) = str2double(gridPreFiles{ii}(indUnd(end-1)+1:indUnd(end)-1));
        if isnan(datePre(ii,1))
            datePre(ii,1) = -9999;
        end
        datePre(ii,2) = str2double(gridPreFiles{ii}(indUnd(end  )+1:indExt(end)-1));
    end
    
elseif regexpbl(fileTyp,'cdf')
    %LOAD ALL NETCDF data:
    %Write required metadata:
    sMeta.mnths = (1:12); %Selects months
    sMeta.crd = crd; %Select bounding box
    sMeta.hisTS = '';
    sMeta.yrsOut = yrs;
    sMeta.yrsClim = [NaN, NaN]; 
    
    %Load all months of data:
    for ii = 1 : 12
        sMeta.currTime = [NaN,sMeta.mnths(ii)];
        sMeta.currVar = 'tmp';
        sDataTmpTemp = read_geodata(gridTmpFiles, sMeta, 'ts');
        sMeta.currVar = 'pre';
        sDataPreTemp = read_geodata(gridPreFiles, sMeta, 'ts');
        
        if ii == 1
            sDataTmp = sDataTmpTemp;
            sDataPre = sDataPreTemp;
        else
            sDataTmp.data = cat(1,sDataTmp.data,sDataTmpTemp.data);
            sDataTmp.time = cat(1,sDataTmp.time,sDataTmpTemp.time);
            sDataPre.data = cat(1,sDataPre.data,sDataPreTemp.data);
            sDataPre.time = cat(1,sDataPre.time,sDataPreTemp.time);
        end 
    end
    
    %Sort:
    [sDataTmp.time, indTmp] = sort(sDataTmp.time,'ascend');
    [sDataPre.time, indPre] = sort(sDataPre.time,'ascend');
    sDataTmp.data = sDataTmp.data(indTmp,:,:);
    sDataPre.data = sDataPre.data(indPre,:,:);
    
    %Find month-yr corresponding to sData.time:
    strCalPre =  NC_cal(sDataPre.attTime);
    [gcmRefPre, gcmUnitsPre] = NC_time_units(sDataPre.attTime);
    if isempty(regexpi(gcmUnitsPre,'days since'))
        error('snowfall_calc:refDate',[gcmUnitsPre ' is an unknown reference system.']);
    end
    datePre = days_2_date(sDataPre.time, gcmRefPre, strCalPre);
    datePre = datePre(:,1:2);
    
    strCalTmp =  NC_cal(sDataTmp.attTime);
    [gcmRefTmp, gcmUnitsTmp] = NC_time_units(sDataTmp.attTime);
    if isempty(regexpi(gcmUnitsTmp,'days since'))
        error('snowfall_calc:refDate',[gcmUnitsTmp ' is an unknown reference system.']);
    end
    dateTmp = days_2_date(sDataTmp.time, gcmRefTmp, strCalTmp);
    dateTmp = dateTmp(:,1:2);
else
    error('snowfall_calc:fileType',[fileTyp ' is an unknown file type, which has not been coded for.'])
end


%Find correspondance between tmp and pre data:
[~, interTmp, interPre] = intersect(dateTmp,datePre,'rows');
vecInterTmp = zeros(length(dateTmp(:,1)), 1);
    vecInterTmp(interTmp) = 1;
vecInterPre = zeros(length(datePre(:,1)), 1);
    vecInterPre(interPre) = 1;

%Remove time-series elements not common to both:
dateTmp(~vecInterTmp,:) = [];
datePre(~vecInterPre,:) = [];

%Find chronological ordering:
[~, indTmp] = sortrows(dateTmp,[1 2]);
[~, indPre] = sortrows(datePre,[1 2]);

%Remove grid elements not common to both meteorological datasets:
if regexpbl(fileTyp,'asc')
    gridTmpFiles(~vecInterTmp) = [];
    gridPreFiles(~vecInterPre) = [];
    gridTmpFiles = gridTmpFiles(indTmp);
    gridPreFiles = gridPreFiles(indPre);
 
    nFiles = numel(gridTmpFiles);
%     foldExt = gridPreFiles{1};
%     indUnd = regexpi(foldExt,'_');
%     foldExt = strrep(foldExt(1:indUnd(end-1)-1),'pre_','snow_');
%     foldExt = strrep(foldExt,'pr_','snow_');

    %Add bias-correction tag or show warning if one bias-corrected and
    %other not:
    if ~isequal(~isempty(regexpi(foldPre,'BC')),~isempty(regexpi(foldTmp,'BC')))
        warning('snow_calc:biasCorrection',['One of the input ' ...
            'datasets appears to be bias-corrected and the other ' ...
            'appears to not be bias-corrected.']);
    end
elseif ~isempty(regexpi(fileTyp,'cdf'))
    sDataTmp.data(~vecInterTmp,:,:) = [];
    sDataPre.data(~vecInterPre,:,:) = [];
    
    nFiles = numel(sDataPre.time);
%     [~,foldExt,~] = fileparts(filePre);
%     indUnd = regexpi(foldExt,'_');
%     foldExt = strrep(foldExt(1:indUnd(end-1)-1),'pr_','');
%     foldExt = strrep(foldExt,'anom_','');
end


%If no files located, return to main function
if nFiles == 0
    warning('snowfall:noFiles',['Snowfall will not be modeled for any '...
        'files.  Two likely causes are (1) an error in one of the '...
        'input folder names or (2) no files in the specified input folder']);
    return
end

%

%%Create Output directory name for time-series files:
indSep = regexpi(foldPre,filesep);
foldSnow = strrep(foldPre(indSep(end-1)+1:end),'_pre_','_snow_');
    foldSnow = strrep(foldSnow,'_pr_','_snow_');
pathSnow = [foldPre(1:indSep(end-1)), foldSnow];
%Create output directory
if ~exist(pathSnow, 'dir')
    mkdir(pathSnow);
end

%Create stats header
if regexpbl(snowMod,'step')
    cellHdr = {'Year','Month', ['# Cells (<=' num2str(param) ' C)'], ['# Cells (>' num2str(param) ' C)']...
        'Avg Temp (entire; deg C)', ['Avg Temp (<=' num2str(param) ' C)'], ['Avg Temp (>' num2str(param) ' C)'], ...
        'Avg Pre (entire; mm)', ['Avg Snow (<=' num2str(param) ' mm)'], ['Avg Snow (>' num2str(param) ' mm)'], 'Tot Snow (m^3)'};
elseif regexpbl(snowMod,'ramp')
    cellHdr = {'Year','Month', ['# Cells (<=' num2str(param(1)) ' C)'], ['# Cells (>' num2str(param(2)) ' C)']...
        'Avg Temp (entire; deg C)', ['Avg Temp (<=' num2str(param(1)) ' C)'], ['Avg Temp (>' num2str(param(2)) ' C)'], ...
        'Avg Pre (entire; mm)', ['Avg Snow (<=' num2str(param(1)) ' mm)'], ['Avg Snow (>' num2str(param(2)) ' mm)'], 'Tot Snow (m^3)'};
end
strHdr = blanks(0);

for ii = 1 : numel(cellHdr)
    if ii ~= numel(cellHdr)
        strHdr = [strHdr char(cellHdr{ii}) ', '];
    else
        strHdr = [strHdr char(cellHdr{ii}) char(10)];
    end
end

%Initialize output matrix:
statsOut = nan(nFiles,numel(cellHdr));

%Write dates to output stats array:
statsOut(:,1:2) = dateTmp(indTmp,1:2);

for ii = 1 : nFiles
    if ~isempty(regexpi(fileTyp,'asc'))
        %Read each file in chronological order::
        [dataCurrTmp, hdrTmp, metaTmp] = read_ESRI(fullfile(foldTmp, gridTmpFiles{ii}));
        [dataCurrPre, hdrPre, metaPre] = read_ESRI(fullfile(foldPre, gridPreFiles{ii}));
        
        [latTmp, lonTmp] = ESRI_hdr2geo(hdrTmp, metaTmp);
        [latPre, lonPre] = ESRI_hdr2geo(hdrPre, metaPre);
        
        rootNm = gridPreFiles{ii};
    elseif ~isempty(regexpi(fileTyp,'cdf'))
        dataCurrTmp = squeeze(sDataTmp.data(ii,:,:));
        dataCurrPre = squeeze(sDataPre.data(ii,:,:));
        latTmp = sDataTmp.lat;
        lonTmp = sDataTmp.lon;
        latPre = sDataPre.lat;
        lonPre = sDataPre.lon;
        
        [~,nmTemp,~] = fileparts(filePre);
        indU = regexpi(nmTemp,'_');
        nmTemp = nmTemp(1:indU(end)-1);
        rootNm = [nmTemp '_' num2str(datePre(ii,1)) '_' num2str(datePre(ii,2)) '.asc'];
        
        hdrTmp = ESRI_hdr(lonTmp, latTmp, 'cor');
    end
    
    %Display progress:
    warning('off','MATLAB:gui:latexsup:UnableToInterpretTeXString');
    fracComplt = round2(100*ii/nFiles,1);
    disp(['Snowfall processing is ' num2str(fracComplt) '% complete.  Data from '...
        char(39) rootNm char(39) ' is being read and processed.'])
    warning('on','MATLAB:gui:latexsup:UnableToInterpretTeXString');
    
    currSnow = zeros(size(dataCurrTmp));
    
    if ~isequal(latTmp,latPre) || ~isequal(lonTmp,lonPre)
        error('snowfall_calc:gridAlign',['The geographical coordiantes for '...
            'precip and temp do not align.']);
    end
       
    %Calculate area (used for measuring snowfall):
    area = area_geodata(lonTmp, latTmp);
    
    if regexpbl(snowMod,'step')
        %Find cells at or below t_thresh Deg C
        indFrozen = find(dataCurrTmp <= param(1));
        indThaw = (1:1:numel(dataCurrTmp))';
        indThaw(indFrozen) = [];
    
        %Number Cells (<= t_thresh C)
        statsOut(ii,3) = sum(indFrozen);

        %Number Cells (> t_thresh C)
        statsOut(ii,4) = numel(indThaw);

        %Avg Temp (entire)
        statsOut(ii,5) = nansum(nansum(area.*dataCurrTmp)) / nansum(nansum(area));

        %Avg Temp (<= t_thresh C)
        if ~isempty(indFrozen)
            statsOut(ii,6) = nansum(nansum(area(indFrozen).*dataCurrTmp(indFrozen))) / nansum(nansum(area(indFrozen)));
        else
            statsOut(ii,6) = 0;
        end

        %Avg Temp (> t_thresh C)
        if ~isempty(indThaw)
            statsOut(ii,7) = nansum(nansum(area(indThaw).*dataCurrTmp(indThaw))) / nansum(nansum(area(indThaw)));
        else
            statsOut(ii,7) = 0;
        end

        %Avg Pre (entire)
        statsOut(ii,8) = nansum(nansum(area.*dataCurrPre)) / nansum(nansum(area));

        %Avg Snow (<= t_thresh C)
        if ~isempty(indFrozen)
            statsOut(ii,9) = nansum(nansum(area(indFrozen).*dataCurrPre(indFrozen))) / nansum(nansum(area(indFrozen)));
        else
            statsOut(ii,9) = 0;
        end

        %Avg Pre (> t_thresh C)
        if ~isempty(indThaw)
            statsOut(ii,10) = nansum(nansum(area(indThaw).*dataCurrPre(indThaw))) / nansum(nansum(area(indThaw)));
        else
            statsOut(ii,10) = 0;
        end

        %Total Snow (m^3)
        if ~isempty(indFrozen)
            statsOut(ii,11) = nansum(nansum(area(indFrozen).*dataCurrPre(indFrozen))) / 1000; %Pre is in units of mm
        else
            statsOut(ii,11) = 0;
        end

        %Calculate Snow
        if ~isempty(indFrozen)
            currSnow(indFrozen) = dataCurrPre(indFrozen); %units are mm (depth)
        else
            currSnow(indFrozen) = 0;
        end
    elseif regexpbl(snowMod,'ramp')
        %Find cells at or below 0 Deg C
        indFrozen = find(dataCurrTmp <= param(1));
        indPart = find(dataCurrTmp > param(1) & dataCurrTmp <= param(2));
        indThaw = (1:1:numel(dataCurrTmp))';
        indThaw([indFrozen;indPart]) = [];
    
        %Calculate Snow
        if ~isempty(indFrozen)
            currSnow(indFrozen) = dataCurrPre(indFrozen); %units are mm (depth)
        end
        if ~isempty(indPart)
            currSnow(indPart) = (param(2) - dataCurrTmp(indPart))/(param(2)-param(1)) ...
                .* dataCurrPre(indPart);
        end
        
        
        %Number Cells (<= t_s C)
        statsOut(ii,3) = sum(indFrozen);

        %Number Cells (> t_r C)
        statsOut(ii,4) = numel(indThaw);

        %Avg Temp (entire)
        statsOut(ii,5) = nansum(nansum(area.*dataCurrTmp)) / nansum(nansum(area));

        %Avg Temp (<= t_s C)
        if ~isempty(indFrozen)
            statsOut(ii,6) = nansum(nansum(area(indFrozen).*dataCurrTmp(indFrozen))) / nansum(nansum(area(indFrozen)));
        else
            statsOut(ii,6) = 0;
        end

        %Avg Temp (> t_r C)
        if ~isempty(indThaw)
            statsOut(ii,7) = nansum(nansum(area(indThaw).*dataCurrTmp(indThaw))) / nansum(nansum(area(indThaw)));
        else
            statsOut(ii,7) = 0;
        end

        %Avg Pre (entire)
        statsOut(ii,8) = nansum(nansum(area.*dataCurrPre)) / nansum(nansum(area));

        %Avg Snow (<= t_s C)
        if ~isempty(indFrozen) || ~isempty(indPart) 
            statsOut(ii,9) = nansum(nansum(area([indFrozen;indPart]).*currSnow([indFrozen;indPart]))) / nansum(nansum(area([indFrozen;indPart])));
        else
            statsOut(ii,9) = 0;
        end

        %Avg Pre (> t_r C)
        if ~isempty(indThaw)
            statsOut(ii,10) = nansum(nansum(area(indThaw).*dataCurrPre(indThaw))) / nansum(nansum(area(indThaw)));
        else
            statsOut(ii,10) = 0;
        end

        %Total Snow (m^3)
        if ~isempty(indFrozen) || ~isempty(indPart) 
            statsOut(ii,11) = nansum(nansum(area([indFrozen;indPart]).*currSnow([indFrozen;indPart]))) / 1000; %Pre is in units of mm
        else
            statsOut(ii,11) = 0;
        end
    end
    
    %Write snowfall data to file:
	if wrtTs == 1 
        fileSnow = strrep(rootNm,'_pre_','_snow_');
        fileSnow = strrep(fileSnow,'_pr_','_snow_');
        fileSnow = strrep(fileSnow,'_anom','');
        write_ESRI_v4(currSnow, hdrTmp, fullfile(pathSnow, fileSnow), 0);
	end
end


%%WRITE STATISTICS FILE:
%Create output filename:
indUnd = regexpi(rootNm,'_');
fileStat = [strrep(rootNm(1:indUnd(end-1)-1),'_pre',''), '_stats.txt'];
fileStat = strrep(fileStat,'pre_','');
fileStat = strrep(fileStat,'pr_','');
fileStat = strrep(fileStat,'_pr','');

pathStatsOut = fullfile(pathSnow, fileStat);
disp(['The spatially aggregated output snow time-series file will be written to ' char(39) ...
    pathStatsOut char(39) '.' char(10) 'The spatially aggregated units '...
    'of the snow are m^3 of water equivalent, while the units for '...
    'spatially-distributed snow are mm of water equivalent.' char(10)]);

%Write Header:
fileID = fopen(pathStatsOut,'w');
fwrite(fileID, strHdr, 'char');
%Write data:
dlmwrite(pathStatsOut,statsOut,'-append');
fclose(fileID);