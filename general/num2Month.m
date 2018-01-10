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

function str = num2Month(n,varargin)

str = cell(numel(n),1);

for ii = 1 : numel(n)
    switch n(ii)
        case 1
            str{ii} = 'January';
        case 2
            str{ii} = 'February';
        case 3
            str{ii} = 'March';
        case 4
            str{ii} = 'April';
        case 5
            str{ii} = 'May';
        case 6
            str{ii} = 'June';
        case 7
            str{ii} = 'July';
        case 8
            str{ii} = 'August';
        case 9
            str{ii} = 'September';
        case 10
            str{ii} = 'October';
        case 11
            str{ii} = 'November';
        case 12
            str{ii} = 'December';
        otherwise
            str{ii} = 'Unknown';
    end

    if ~isempty(varargin)
       str{ii} = str{ii}(1:varargin{1}); 
    end
end    