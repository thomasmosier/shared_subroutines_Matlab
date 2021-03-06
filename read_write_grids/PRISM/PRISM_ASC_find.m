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

function [filesUse, varargout] = PRISM_ASC_find(dirTS, varargin)



%Designed to work with ASCII data files for which the date is indicated in
%the file as '*year_mnth.asc'



if numel(varargin) == 2
   yrLoad = varargin{1};
   mnthLoad = varargin{2}; 
elseif numel(varargin) == 1
    yrLoad(1) = min(varargin{1}(:,1));
    yrLoad(2) = max(varargin{1}(:,1));
    mnthLoad = varargin{1}(:,2);
else
    yrLoad = nan(1,2);
    mnthLoad = (1:12);
end
mnthsUniq = unique(mnthLoad);
nMnths = numel(mnthsUniq);


if all(~isnan(mnthLoad)) && all(~isnan(yrLoad)) 
    if numel(mnthLoad) ~= numel(yrLoad) && nMnths == numel(mnthLoad)  
        if numel(yrLoad) == 2
            nYrsLd  = yrLoad(2) - yrLoad(1) + 1;
        else
            nYrsLd = numel(unique(yrLoad));
        end

        ldTs = nan(nYrsLd*nMnths,2);

        for ii = 1 : nMnths
            ldTs((ii-1)*nYrsLd+1 : ii*nYrsLd,1) = (yrLoad(1):yrLoad(end))';
            ldTs((ii-1)*nYrsLd+1 : ii*nYrsLd,2) = mnthsUniq(ii);
        end
    elseif numel(yrLoad) == numel(mnthLoad)
        ldTs(:,1) = yrLoad;
        ldTs(:,2) = mnthLoad;
    else
        error('ASC_file_find:dateFormat','The input data format is unknown.')
    end
else
   ldTs = nan(1,2);  
end


%Find root filename and date range:
fileTemp = dir([dirTS filesep '*.asc']);
filesTS = cell(length(fileTemp),1); 
for ii = 1 : length(fileTemp)
    filesTS(ii) = cellstr(fileTemp(ii).name);
end

%Find dates of all files:
ind4Digit = regexpi(filesTS{1},'\d{4}');
ind6Digit = regexpi(filesTS{1},'\d{6}');

if isempty(ind4Digit) && isempty(ind6Digit)
    error('PrismAscFind:unknownType',['There are no numbers with 4 or 6 digits in the file ' filesTS{1} '. The file type cannot be determined.'])
elseif ~isempty(ind4Digit) && ~isempty(ind6Digit)
    if ind4Digit(end) > ind6Digit(end)
        type = 'clim';
    else
        type = 'ts';
    end
elseif isempty(ind4Digit) && ~isempty(ind6Digit)
    type = 'ts';
else
    type = 'clim';
end


if strcmpi(type, 'ts')
    fileDate = nan(length(filesTS), 2);
    for jj = length(filesTS) : -1 : 1
        ind6Digit = regexpi(filesTS{jj},'\d{6}');

        if isempty(ind6Digit)
            filesTS(jj) = [];
            fileDate(jj,:) = [];
        else
            fileDate(jj,:) = [str2double(filesTS{jj}(ind6Digit:ind6Digit+3)), str2double(filesTS{jj}(ind6Digit+4:ind6Digit+5))];
        end
    end
else
    error('PrismAscFind:programClim','Climatological file formats have not been programmed for.')
end


%Find dates requested:
if all(~isnan(ldTs))
    [~, indUse, ~] = intersect(fileDate,ldTs,'rows');
else
    indUse = (1:numel(fileDate(:,1)));
end

%Find file names of all files to use:
filesUse = cell(numel(indUse),1);
for ii = 1 : numel(indUse)
    filesUse(ii) = fullfile(dirTS, cellstr(filesTS(indUse(ii))));
end
fileDate = fileDate(indUse,:);
   
if any2d(~isnan(ldTs)) && numel(indUse) ~= sum(~isnan(ldTs(:,1)))
    warning('prismAscFind:nWrong',[num2str(numel(indUse)) ...
        ' have been loaded from the folder ' char(39) strrep(dirTS,filesep,':') char(39) ...
        ' but the reqested number to load is ' num2str(numel(ldTs(:,1))) '.']);
end

%Put data in order:
[fileDate, indSort] = sortrows(fileDate);
filesUse(indSort);


if numel(varargin) == 1
    [fileDate, ~, indUse] = intersect(varargin{1},fileDate,'rows');
    filesUse = filesUse(indUse);
end

if nargout > 1
   varargout{1} = fileDate; 
   if nargout > 2
       if numel(indUse) > 0
            [~, varargout{2}, metaESRI] = read_ESRI(filesUse{1});
       else
           varargout{2} = nan(6,1);
       end
       
        if nargout > 3
            if numel(indUse) > 0
                varargout{3} = metaESRI;
            else
                varargout{3} = cell(6,1);
            end
        end
   end
end
