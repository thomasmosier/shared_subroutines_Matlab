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

function dateFilled = date_vec_fill(dateStart,dateEnd,cal)
%Assumes Gregorian calendar

if ~regexpbl(cal,{'greg','365_day','no_leap'})
   error('date_vec_fill:unkownCal',[cal ' has not been programmed for.  '...
       'Currently this function only works with the Gregorian Calendar.']); 
end


if numel(dateStart(1,:)) ~= numel(dateEnd(1,:))
    if numel(dateStart(1,:)) < numel(dateEnd(1,:))
        strCrop = 'starting';
        dateStart = dateStart(1:numel(dateEnd(1,:))); 
    else
        strCrop = 'ending';
        dateEnd = dateEnd(1:numel(dateStart(1,:))); 
    end

    warning('date_vec_fill:diffLength',['The time resolution of the '...
        'input vectors is different.  The vectors are being cropped '...
        'to the resolution of the ' strCrop ' vector.']);
end

if numel(dateStart(1,:)) == 1
    nSteps = dateEnd - dateStart + 1;
elseif numel(dateStart(1,:)) == 2
    nSteps = 12*(dateEnd(1) - dateStart(1) + 1);
elseif numel(dateStart(1,:)) == 3
    nSteps = 366*(dateEnd(1) - dateStart(1) + 1);
elseif numel(dateStart(1,:)) == 4
    nSteps = 8784*(dateEnd(1) - dateStart(1) + 1);
end


dateFilled = nan(nSteps, numel(dateStart));
dateFilled(1,:) = dateStart;

if numel(dateStart(1,:)) == 1
    for ii = 1 : dateEnd - dateStart 
        dateFilled(ii+1,:) = dateStart + ii;
    end
elseif numel(dateStart(1,:)) == 2
    ii = 1;
    while ~isequal(dateFilled(ii,:), dateEnd)
        if dateFilled(ii,2) ~= 12
           dateFilled(ii+1,:) = dateFilled(ii,:) + [0,1];
        else
            dateFilled(ii+1,:) = [dateFilled(ii,1) + 1, 1];
        end

        ii = ii + 1;
    end
elseif numel(dateStart(1,:)) == 3
    ii = 1;
    while ~isequal(dateFilled(ii,:), dateEnd)
        if dateFilled(ii,2) == 2
            if regexpbl(cal, 'gregorian') && is_leap_v2(dateFilled(ii,1))
                dayTest = 29;
            else
                dayTest = 28;
            end
        else
            dayTest = eomday(dateFilled(ii,1), dateFilled(ii,2));
        end
        
        if dateFilled(ii,3) ~= dayTest
           dateFilled(ii+1,:) = dateFilled(ii,:) + [0,0,1];
        else
            if dateFilled(ii,2) == 12
                dateFilled(ii+1,:) = [dateFilled(ii,1) + 1, 1, 1];
            else
                dateFilled(ii+1,:) = [dateFilled(ii,1), dateFilled(ii,2) + 1, 1];
            end
        end

        ii = ii + 1;
    end
elseif numel(dateStart(1,:)) == 4
    ii = 1;
    while ~isequal(dateFilled(ii,:), dateEnd)
        if dateFilled(ii,4) ~= 24;
            dateFilled(ii+1,:) = dateFilled(ii,:) + [0,0,0,1];
        else
            if dateFilled(ii,3) ~= eomday(dateFilled(ii,1),dateFilled(ii,2))
               dateFilled(ii+1,:) = [dateFilled(ii,1:3) + [0,0,1], 1];
            else
                if dateFilled(ii,2) == 12
                    dateFilled(ii+1,:) = [dateFilled(ii,1) + 1, 1, 1, 1];
                else
                    dateFilled(ii+1,:) = [dateFilled(ii,1), dateFilled(ii,2) + 1, 1, 1];
                end
            end
        end

        ii = ii + 1;
    end
end

dateFilled( isnan(dateFilled(:,1)), :) = [];