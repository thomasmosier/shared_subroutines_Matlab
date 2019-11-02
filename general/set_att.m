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

function listAtt = set_att(listAtt, strAtt, val)
    
ind = strcmpi(listAtt,strAtt);

if sum2d(ind) == 1
    [row, col] = find(ind == 1);
    listAtt{row, col+1} = val;
elseif sum2d(ind) > 1
    error('set_att:mult', ['Multiple attributes labeled ' char(39) strAtt char(39) ' were found.']);
elseif sum2d(ind) == 0
    warning('set_att:none', ['No attributes labeled ' char(39) strAtt char(39) ' were found.']);
else
    warning('set_att:none', 'Unknown value returned from strcmpi.');
end