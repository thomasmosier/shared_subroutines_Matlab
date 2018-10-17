function [dataMat, dates, hdrOut, metaHdr, bndsDdAph, bndsIndAph] = read_APHRODITE_bin(pathAphro, gridDd, mnthsLd, cropSide)


%APHRODITE data stored as double
type = 'float';


% hdrEsri = [ ncol; ...     %ncols
%             nrow; ...     %nrows
%             llx; ...  %xllcorner
%             lly; ...    %yllcorner
%             xStep; ...  %cellsize
%             noData];    %NODATA_value


%Load data:
fidAphData = fopen(pathAphro); %Open binary APHRODITE file. 
aphVec = fread(fidAphData, Inf, type, 0, 'l');  %Read APHRODITE binary file.
fclose(fidAphData);
nData = numel(aphVec);


noDataOut = -9999;
%Determine region and coordinate properties:
if regexpbl(pathAphro, '_MA_')
    %Lon 60.0E - 150.0E, Lat 15.0S - 55.0N
    hdrRaw = [ nan; ...     %ncols
        nan; ...     %nrows
        60; ...  %xllcorner
        -15; ...    %yllcorner
        nan; ...  %cellsize
        noDataOut];    %NODATA_value
    
    %0.25 res: 360*280;
    nLat25 = 280;
    nLon25 = 360;
    %0.5 res: 180*140;
    nLat5 = 140;
    nLon5 = 180;
    
    n25 = nLat25*nLon25; 
    n5   = nLat5* nLon5;

    %Determine resolution (either 0.25 or 0.5)
    if 0.5*nData / n25 > 364 && 0.5*nData / n25 < 367
        hdrRaw(1) = nLon25; %ncols
        hdrRaw(2) = nLat25; %nrows
        hdrRaw(5) = 0.25; %stepsize
    elseif 0.5*nData / n5 > 364 && 0.5*nData / n5 < 367
        hdrRaw(1) = nLon5; %ncols
        hdrRaw(2) = nLat5; %nrows
        hdrRaw(5) = 0.5; %stepsize
    else
        error('readAphroditeData:unknownResolution', ['The current APHRODITE ' ...
            'data file for Monsoon Asia does not appear to have a resolution of 0.25 or 0.5 '...
            'decimal degrees.']);
    end
elseif regexpbl(pathAphro, '_JP_')
    %Lon 123.0E-146.0E; Lat 24.0N-46.0N
    hdrRaw = [ nan; ...     %ncols
        nan; ...     %nrows
        123; ...  %xllcorner
        24; ...    %yllcorner
        nan; ...  %cellsize
        noDataOut];    %NODATA_value
    
    %0.05 res: 460*440
    nLat05 = 440;
    nLon05 = 460;

    n05   = nLat05* nLon05;
    %Determine resolution (either 0.25 or 0.5)
    if 0.5*nData / n05 > 364 && 0.5*nData / n05 < 367
        hdrRaw(1) = nLon05; %ncols
        hdrRaw(2) = nLat05; %nrows
        hdrRaw(5) = 0.05; %stepsize
    else
        error('readAphroditeData:unknownResolution', ['The current APHRODITE ' ...
            'data file for Japan does not appear to have a resolution of 0.05 '...
            'decimal degrees.']);
    end
else
    error('readAphroditeData:unknownRegion', ['The region for file ' pathAphro ' has not been programmed.'])
end


% Each element (both temperature and
% rain gauge information) is written as a 4-byte floating-point number
% in little endian byte order.  Users should swap the byte order to
% big endian if necessary.  There are no 'space', 'end of record', or
% 'end of file' marks in between.  In the case of the 0.5-degree APHRO_TAVE_MA
% product, the size of a file (0.5-degree grid) is
%    4 bytes x 180 x 140 x 2 fields x 365 days = 73,584,000 bytes
% for a non-leap year, or 73,785,600 bytes for a leap year.
    

%Reshape input based on above determinations:
dataMat = single(permute(reshape(aphVec,  hdrRaw(1), hdrRaw(2), 2, []), [4,2,1,3]));
%Remove station contributions:
dataMat = squeeze(dataMat(:,:,:,1));
%Flip latitudinal orientation
dataMat = flip(dataMat,2);

% figure;
% imagesc(squeeze(dataMat(100,:,:))); colorbar;

% figure;
% imagesc(squeeze(dataMat(200,:,:))); colorbar;

%Crop lat/lon
if ~isempty(gridDd)
    [bndsDdAph, bndsIndAph] = adj_bounds(hdrRaw, gridDd, cropSide);    %Find indices to crop APHRODITE data to.
    dataMat = dataMat(:, bndsIndAph(4):bndsIndAph(3), bndsIndAph(1):bndsIndAph(2));  %Crop APHRODITE to desired region.
else
    bndsDdAph = [hdrRaw(3) + hdrRaw(5)/2, hdrRaw(3) + hdrRaw(5)*(0.5 + numel(dataMat(1,1,:))), hdrRaw(4) - hdrRaw(5)/2, hdrRaw(4) - (numel(dataMat(1,:,1)) + 0.5)*hdrRaw(5)];
    bndsIndAph = [1, numel(dataMat(1,1,:)), 1, numel(dataMat(1,:,1))];
end
  

%Set no data to NaN
dataMat(dataMat < -80) = nan;
 
% figure;
% imagesc(squeeze(dataMat(200,:,:))); colorbar;

%Find year
yrIn = str2double(pathAphro(end-3:end));

%Create date vector for full year
dates = nan(366,3, 'single');
dates(:,1) = yrIn;
cntr = 0;
for ii = 1 : 12
    daysCurr = eomday(yrIn, ii);
    
    indCurr = (cntr + 1 : cntr + daysCurr);
    
    dates(indCurr,2) = ii;
    dates(indCurr,3) = (1:daysCurr);
    
    cntr = cntr + daysCurr;
end

%Remove last row if not leap year
if cntr == 365
    dates(366,:) = [];
end


%Remove dates and data not in months load
monthsRem = setdiff((1:12), mnthsLd);
indDateRem = ismember(dates(:,2), monthsRem);
if any(indDateRem)
    dates(indDateRem,:) = [];
    dataMat(indDateRem,:,:) = [];
end


%Write output header
hdrOut(1) = length(dataMat(1,1,:));
hdrOut(2) = length(dataMat(1,:,1));
if ~isempty(gridDd)
    hdrOut(3) = bndsDdAph(1) - hdrRaw(5)/2;
    hdrOut(4) = bndsDdAph(3) - hdrRaw(5)/2;
else
    hdrOut(3) = hdrRaw(3);
    hdrOut(4) = hdrRaw(4);
end
hdrOut(5) = hdrRaw(5);
hdrOut(6) = hdrRaw(6);

metaHdr = {'NCOLS';'NROWS';'XLLCORNER';'YLLCORNER';'CELLSIZE';'NODATA_VALUE'};