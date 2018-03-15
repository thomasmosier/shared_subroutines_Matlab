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

function cal =  NC_cal(attTime)

strNone = 'unknown';
if isempty(attTime)
    cal = strNone;
else
    tCalInd = strcmpi(attTime,'calendar');
    
    if all(tCalInd == 0)
        cal = strNone;
    else
        [tCalRow, tCalCol] = find(tCalInd == 1);
        cal = attTime{tCalRow, tCalCol+1};
    end
end

if strcmpi(cal, strNone)
    warning('ncCal:unknownCal',['The calendar for the current file is '...
        'not known. It is therefore being assumed to be Gregorian.']);
    cal = 'Gregorian';
end