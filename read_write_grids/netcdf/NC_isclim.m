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

function [boolClim, yrs, mnth] = NC_isclim(sData)

if isfield(sData, 'attTime')
    indType = strcmpi(sData.attTime,'type');
else
    indType = 0;
end

if sum2d(indType) == 0
    boolClim = 0;
else
    [tTypRow, tTypCol] = find(indType == 1);
    strType = sData.attTime{tTypRow,tTypCol+1};

    boolClim = ~isempty(regexpi(strType,'clim'));
end

if boolClim == 1
    numUse = str2double(regexpi(strType,'\d*','match'));
    yrsUse = numUse(numUse > 999);
    yrs = sort(yrsUse);
    
    mnthUse = numUse(numUse < 13);
    if length(mnthUse) == 2
        if mnthUse(1) == mnthUse(2)
            mnth = mnthUse(1);
        else
            warning('NC_iscell:inconsistentMonths',['The two months '...
                'present in the climatology string are ' ...
                num2str(mnthUse(1)) ' and ' num2str(mnthUse(2)) '. ' ...
                num2str(mnthUse(1)) ' is being used.']);
            mnth = mnthUse(1);
        end
    elseif length(mnthUse) == 1
        mnth = mnthUse;
    end
else
    yrs = nan(1,2);
    mnth = nan;
end