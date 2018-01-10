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


function files = find_files(fold,ext,varargin) 
%Find files in fold, with 'ext'.  'varargin' can be string to return only a
%subset of the files (i.e. if the file names contain the string)


filesTemp = dir(fullfile(fold, ['*.' ext]));
    files = struct2cell(filesTemp);
    files = files(1,:);
    
if ~isempty(varargin(:)) && ~isempty(varargin{1})
    for ii = numel(files(:)) : -1 : 1
       if ~regexpbl(files{ii},varargin{1}) 
            files(ii) = [];
       end
    end
end