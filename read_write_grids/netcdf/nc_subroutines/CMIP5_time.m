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

function varargout = CMIP5_time(fileData)



%Find yrs and mnths of CMIP5 GCM data and CRU data based upon name string.


%Initialize output:
yrsGcm = nan(2,1);
mnthsGcm = nan(2,1);
daysGcm = nan(2,1);


if exist(fileData,'file')
    [~, fileData, ~] = fileparts(fileData);
elseif strcmpi(fileData(end-2),'.')
    fileData = fileData(1:end-3);
elseif strcmpi(fileData(end-3),'.')
    fileData = fileData(1:end-4);
end

indDate = regexpi(fileData,'_');
indPer = regexpi(fileData,'\.');
indHyp = regexpi(fileData,'-');

if ~isempty(indDate) && isempty(indPer) && ~isempty(indHyp) %CMIP5 NetCDF format
    strDate = fileData(indDate(end)+1:end);
    indHyp = regexpi(strDate,'-');
    yrsGcm = [str2double(strDate(1:4)); str2double(strDate(indHyp +1 : indHyp +4))];
    mnthsGcm = [str2double(strDate(5:6)); str2double(strDate(indHyp +5 : indHyp +6))];
    if ~isnan(str2double(strDate(7)))
        daysGcm = [str2double(strDate(7:8)); str2double(strDate(indHyp +7 : indHyp +8))];
    end
elseif ~isempty(indDate) && isempty(indPer) && isempty(indHyp) %Maybe ASCII format?
    %Find contiguous numbers that occur after underscores:
    indNumb = regexpi(fileData, '\d');
    indFrst = regexpi(fileData, '\d+');
    dates = regexpi(fileData,'(\d+)','tokens');
    if ~isempty(dates)
        indGap = [1, find(abs(diff(indNumb)) > 1)];
        if numel(indGap) > 1
            cntr = 0;
            for ii = 2 : numel(indGap)
                if indGap(ii- cntr) - indGap(ii-1- cntr) > 5
                    indGap(ii - 1 - cntr) = [];
                    dates(ii - 1 - cntr) = [];
                    cntr = cntr + 1;
                end
            end
        else
            if iscell(dates{1})
                dates{1} = cell2mat(dates{1});
            end
            if str2double(dates{1}) > 32 %Likely this is a year file:
                yrsGcm = [str2double(dates{1}); str2double(dates{1})];
                mnthsGcm = [1; 12];
                dates{1} = {num2str(yrsGcm(1)); num2str(yrsGcm(1))};
                dates{2} = {'1'; '12'};
            else
                dates{1} = yrsGcm;
                dates{2} = mnthsGcm;
            end
        end

        for ii = numel(dates) : -1 : 1
            if str2double(dates{ii}) > 10^4
                dates{ii} = [];
            end
        end

        if str2double(dates{1}) > 31 %Year first
            yrsGcm(1:2) = str2double(dates{1});
            mnthsGcm(1:2) = str2double(dates{2});
            if numel(dates) > 2
                daysGcm(1:2) = str2double(dates{3});
            end
        elseif str2double(dates{1}) > 12 %Day first
            daysGcm(1:2) = str2double(dates{1});
            if str2double(dates{2}) > 31 %Year 2nd
                yrsGcm(1:2) = str2double(dates{2});

                if numel(dates) > 2
                    mnthsGcm(1:2) = str2double(dates{3});
                end
            else %Day 2nd
                mnthsGcm(1:2) = str2double(dates{2});

                if numel(dates) > 2
                    yrsGcm(1:2) = str2double(dates{3});
                end
            end
        else %Month first
            mnthsGcm(1:2) = str2double(dates{1});

            if str2double(dates{2}) > 31 %Year 2nd
                yrsGcm(1:2) = str2double(dates{2});

                if numel(dates) > 2
                    daysGcm(1:2) = str2double(dates{3});
                end
            else %Day 2nd
                daysGcm(1:2) = str2double(dates{2});

                if numel(dates) > 2
                    yrsGcm(1:2) = str2double(dates{3});
                end
            end

        end
    end
elseif ~isempty(indPer) %Maybe CRU NetCDF data or CFSR?
    indNum = regexpi(fileData,'(\d+)');
    dates = regexpi(fileData,'(\d+)','tokens');
    if ~isempty(dates)
        cntDates = 1;
        mnthsTemp = nan(2,1);
        for ii = 1 : length(dates(:))
            currNum = char(dates{ii});
            if numel(currNum) == 4
                yrsGcm(cntDates) = str2double(currNum);
                cntDates = cntDates + 1;
            elseif numel(currNum) == 6 %|| numel(currNum) == 7
                yrsGcm(cntDates) = str2double(currNum(1:4));
                mnthsGcm(cntDates) = str2double(currNum(5:end));
                cntDates = cntDates + 1;
            elseif numel(currNum) == 8
                yrsGcm(   cntDates) = str2double(currNum(1:4));
                mnthsGcm(cntDates) = str2double(currNum(5:6));
                daysGcm(  cntDates) = str2double(currNum(7:8));
                cntDates = cntDates + 1;
            end
        end
        if cntDates > 3
            yrsGcm = nan(2,1);
        end
        yrsGcm = sort(yrsGcm);
        if sum(~isnan(yrsGcm)) == 2 && yrsGcm(1) ~= yrsGcm(2) && isequal(mnthsTemp,[1;12])
            mnthsGcm = [1;12];
        end
    end
end

%If ending date nan, set to be same as beginning date:
if ~isnan(yrsGcm(1)) && isnan(yrsGcm(2))
    yrsGcm(2) = yrsGcm(1);
end
if ~isnan(mnthsGcm(1)) && isnan(mnthsGcm(2))
    mnthsGcm(2) = mnthsGcm(1);
end
if ~isnan(daysGcm(1)) && isnan(daysGcm(2))
    daysGcm(2) = daysGcm(1);
end

%write special code for GCM with units of 'BP' (before present)
if all2d(isnan(yrsGcm)) && regexpbl(fileData,'BP')
    yrsGcm = [-9999;-9999];
    mnthsGcm = [-9999;-9999];
end
if numel(yrsGcm) ~= 2
    warning('CMIP5_time:nYrs',['This function identified ' ...
        num2str(numel(yrsGcm)) ' that may be year references.  ' ...
        'Using normal CMIP5 naming convention, there should be a ' ...
        'starting year and ending year.']);
end

if nargout <= 1
    varargout{1} = [yrsGcm, mnthsGcm, daysGcm];
elseif nargout == 2
   varargout{1} = yrsGcm;
   varargout{2} = mnthsGcm;
elseif nargout == 3
   varargout{1} = yrsGcm;
   varargout{2} = mnthsGcm;
   varargout{3} = daysGcm; 
end
