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

function bl = regexpbl(listFull, listSub, varargin)


%Function matches strings from 'str2' and 'str1'.
%Default is to return 1 if at least one entry is contained in both.
%'AND' conditional requires all strings from 'str2' be present in 'str1'.

%INPUTS:
%'str1': string or cell array of potential candidate strings
%'str2': string or cell array of strings to match with 'str2'
%'varargin': either 'or' or 'and' (default is 'or').
    %With 'or' any string in listSub must be contained in any single string in listFull.
    %With 'and' each string in listSub must be contained a string in listFull

%OUTPUT:
%'bl': 1 or 0 determined by whether condition was met.

if isempty(listFull)
    bl = 0;
    return
end

chr1 = ischar(listFull);
cll1 = iscell(listFull);
chr2 = ischar(listSub);
cll2 = iscell(listSub);

%Varargin can be used to specify logical 'or' or 'and' comparison 
if ~isempty(varargin) && ~isempty(varargin{1})
    typ = varargin{1};
    
    if isempty(regexpi(typ,'and')) && isempty(regexpi(typ,'or'))
        warning('regexpbl:type',['Type is set to ' char(39) typ char(39) ...
            ', which is not a valid option.  It is therefore being set to ' ...
            char(39) 'or' char(39) '.']);
        typ = 'or';  
    end
else
   typ = 'or'; 
end



%HANDLE SPECIAL CASES (this is faster than full version below):
if chr1 && chr2 %Both are character strings
    bl = ~isempty(regexpi(listFull, listSub));
    return
elseif ~chr1 && cll1 && chr2 && strcmpi(typ, 'or') %arg1 is cell array, arg2 is string and selection is 'or'
    bl = 0;
    for ii = 1 : numel(listFull(:))
        if ischar(listFull{ii}) && ~isempty(regexpi(listFull{ii},listSub))
           bl = 1;
           return
        end
    end
    
    return
elseif chr1 && ~chr2 && cll2 && strcmpi(typ, 'or') %arg1 is string, arg2 is cell array and selection is 'or'
    bl = 0;
    for ii = 1 : numel(listSub(:))
        if ischar(listSub{ii}) && ~isempty(regexpi(listFull, listSub{ii}))
           bl = 1;
           return
        end
    end
    
    return
end


%%HANDLE ALL CASES (SLOWER)
if iscell(listFull)
   nIt1 = length(listFull(:));
else
    nIt1 = 1;
    str1Temp = listFull;
    listFull = cell(1,1);
    listFull{1} = str1Temp;
end

if iscell(listSub)
   nIt2 = length(listSub(:));
else
    nIt2 = 1;
    str2Temp = listSub;
    listSub = cell(1,1);
    listSub{1} = str2Temp;
end

%Initialize output:
bl = 0;
blTemp = zeros(length(listSub(:)),1);

%Loop over all combinations of inputs:
for ii = 1 : nIt1
    for jj = 1 : nIt2
        if ~ischar(listFull{ii}) || ~ischar(listSub{jj})
            continue
        elseif ~isempty(regexpi(listFull{ii},listSub{jj})) %|| ~isempty(regexpi(listSub{jj},listFull{ii}))
            blTemp(jj) = 1;
            if strcmpi(typ,'or')
                bl = 1;
                return
            end
        end
    end
end

if strcmpi(typ,'and') && sum(blTemp) >= numel(listSub(:))
    bl = 1; 
end
