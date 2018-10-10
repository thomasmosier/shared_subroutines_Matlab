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

function strMnth = month_num2str(num, varargin)

typ = 'full';
for ii = 1 : numel(varargin(:))
    if regexpbl(varargin{ii}, 'type')
        typ = varargin{ii+1};
    end
end


strMnth = cell(numel(num),1);

for ii = 1 : numel(num)
    switch num(ii)
        case 1
            strMnth{ii} = 'January';
        case 2
            strMnth{ii} = 'February';
        case 3
            strMnth{ii} = 'March';
        case 4
            strMnth{ii} = 'April';
        case 5
            strMnth{ii} = 'May';
        case 6
            strMnth{ii} = 'June';
        case 7
            strMnth{ii} = 'July';
        case 8
            strMnth{ii} = 'August';
        case 9
            strMnth{ii} = 'September';
        case 10
            strMnth{ii} = 'October';
        case 11
            strMnth{ii} = 'November';
        case 12
            strMnth{ii} = 'December';
        otherwise
            strMnth{ii} = 'Unknown';
    end

    
    if strcmpi(typ, 'first')
        strMnth{ii} = strMnth{ii}(1);
    elseif strcmpi(typ, 'three')
        strMnth{ii} = strMnth{ii}(1:3);
    elseif ~strcmpi(typ, 'full')
        error('monthNum2Str:unknownType', ['Option for type ' typ ' has not been programmed for.']);
    end
end
    