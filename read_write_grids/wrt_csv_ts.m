function wrt_csv_ts(pathData, sData, date, varWrt, varargin)

yrsOut = nan(1,2);
if ~isempty(varargin(:))
    for ii = 1 : numel(varargin(:))
       if strcmpi(varargin{ii}, 'yrs') || strcmpi(varargin{ii}, 'years')
           yrsOut = varargin{ii+1};
       end
    end
end

%Create output header and arrays:
switch numel(date(1,:)) 
    case 1
        hdrOut = 'year';
    case 2
        hdrOut = 'year, month';
    case 3
        hdrOut = 'year, month, day';
    case 4
        hdrOut = 'year, month, day, hour';
    otherwise
        error('writeGeodata:unknownDate',['The date array has ' ...
            num2str(numel(date(1,:))) ' precision. This has not been programmed for.']);
end

if ischar(varWrt)
   varWrt = {varWrt}; 
end

outDataArry = date;
for ii = 1 : numel(varWrt)
    hdrOut = [hdrOut ', ' varWrt{ii}];

    if numel(date(:,1)) ~= numel(sData.(varWrt{ii})(:))
        error('writeGeodata:diffSize', ['The date array has length ' ...
            num2str(numel(date(:,1))) ' and the data array has length ' ...
            num2str(numel(sData.(varWrt{ii})(:))) '.']);
    end

    outDataArry = [outDataArry, sData.(varWrt{ii})];
end

%Crop output years:
if all(~isnan(yrsOut))
   indKeep = find(outDataArry(:,1) >= min(yrsOut) & outDataArry(:,1) <= max(yrsOut));
   
   outDataArry = outDataArry(indKeep, :);
   date = date(indKeep, :);
end

%Make directory:
[foldWrt, fileWrt, ext] = fileparts(pathData);

%Estimate time step
if (numel(date(1,:)) == 2 || all(date(:,3) == date(1,3)))
    timeStep = 'monthly';
elseif numel(date(1,:)) == 3 && mode(abs(diff(date(:,3)))) == 1
    timeStep = 'daily';
else
    timeStep = 'unknown';
end

if numel(date(1,:)) == 2
    dayStart = 0;
    dayEnd = 0;
else
    dayStart = date(1,3);
    dayEnd = max(date(:,3));
end
      
%Remove any existing date from file name:
ind8 = regexpi(fileWrt, '\d{8,8}');
if numel(ind8) == 2
    fileWrt(ind8(1):ind8(2)+7) = [];
end
ind6 = regexpi(fileWrt, '\d{6,6}');
if numel(ind6) == 2
    fileWrt(ind6(1):ind6(2)+5) = [];
end

%Add date to file name:
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

    fileWrt = [fileWrt, '_' num2str(date(1,1)), date(1,2) dayStart '-', ...
        num2str(date(end,1)), date(end,2), dayEnd, ext];
elseif regexpbl(timeStep, 'monthly')
    fileWrt = [fileWrt, '_' num2str(date(1,1)), date(1,2) '-', ...
        num2str(date(end,1)), date(end,2), ext];
elseif regexpbl(timeStep, 'unknown')
    fileWrt = [fileWrt '_unknown-dates'];
else
    error('writeGeodata:unknownDates','The data appear to have an unexpected time resolution.')
end

fileWrt = strrep(fileWrt, '__', '_');

        
if ~exist(foldWrt, 'dir')
    mkdir(foldWrt);
end

%If output file already exists, delete
if exist(pathData, 'file')
    delete(pathData)
end

%Check correct extension
if isempty(ext)
    pathData = [pathData '.csv'];
elseif ~isequal(ext,'.csv')
   pathData = fullfile(foldWrt, [fileWrt, '.csv']); 
end

%Write file
fOut = fopen(pathData,'w+');
fprintf(fOut,'%s\n', hdrOut);
fclose(fOut);

dlmwrite(pathData, outDataArry, '-append');