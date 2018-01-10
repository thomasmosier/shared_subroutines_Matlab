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

function gcmDateVec = time_kaBP_2_standard(gcmTime,gcmRef,cal)
%Function only necessary when NetCDF units are 'ka BP' (thousands of years
%before present).  Assumes monthly data and that calendar is 'no leap'

if regexpbl(cal,{'no','leap'},'and')
    gcmDateVec = nan(numel(gcmTime),3);
    gcmDateVec(:,1) = gcmRef(1) + ceil(10^3*gcmTime); 
    gcmDateVec(:,2) = rem((1:numel(gcmTime)),12);
    gcmDateVec(gcmDateVec(:,2) == 0, 2) = 12;
    gcmDateVec(:,3) = 15;
else
   error('time_kaBP_2_standard:calendar',['No case has been written for calendar of type ' cal]);
end


end