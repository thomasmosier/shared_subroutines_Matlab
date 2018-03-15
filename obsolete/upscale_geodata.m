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

function [s1, s2, nRe] = upscale_geodata(s1, s2, varargin)

warning('upscaleGeodata:obsolete','This function is obsolete. Upgrade to upscale_geodata_v3');

varLon = 'longitude';
varLat = 'latitude';
%Find which of the two grids is at a lower spatial resolution and
%interpolate (with area weighting) the higher resolution data to the lower
%resolution grid

%Inputs:
%'s1' = structure array (Downscaling Package format)
%'s2' = structure array (Downscaling Package format)
%'varargin{1}' = if 1, scales s2 to s1 grid; if 2, scales s1 to s2 grid;
    %default scales to higher-spatial resolution grid
%'varargin{2}' = 'sum' or 'nansum' (default is nansum)

%Outputs:
%nRe = number corresponding to dataset that was upscaled ('NaN' indicates
%neither.
%'.dataRe' = upscaled data corresponding to other datasets lat and lon
%'.(varLat)Re' = lat vector corresponding to upscaled data
%'.(varLon)Re' = longitude vector corresponding to upscaled data

%Make all lat and lon vectors column vectors:
[rD, cD] = size(s1.(varLat));
if rD == 1 && cD ~= 1
    s1.(varLat) = s1.(varLat)';
end
[rD, cD] = size(s1.(varLon));
if rD ~= 1 && cD == 1
    s1.(varLon) = s1.(varLon)';
end
[rD, cD] = size(s2.(varLat));
if rD == 1 && cD ~= 1
    s2.(varLat) = s2.(varLat)';
end
[rD, cD] = size(s2.(varLon));
if rD ~= 1 && cD == 1
    s2.(varLon) = s2.(varLon)';
end

d1Lat = diff(s1.(varLat));
d2Lat = diff(s2.(varLat));
d1Lon = diff(s1.(varLon));
d2Lon = diff(s2.(varLon));

d1Avg = mean([abs(d1Lat)', abs(d1Lon)]);
d2Avg = mean([abs(d2Lat)', abs(d2Lon)]);

%Initialize number to indicate which grid was upscaled
nRe = NaN;

type = 'nansum';
if isequal(s1.(varLat),s2.(varLat)) && isequal(s1.(varLon),s2.(varLon)) %Both have common grid
    %No grid upscaled
    return 
elseif ~isempty(varargin(:)) && ~isempty(varargin{1})
    if varargin{1} == 0 %If 0, don't regrid either
        return
    elseif varargin{1} == 1 %If 1, regrid s2 to s1
        nRe = 2;
    elseif varargin{1} == 2 %If 2, regrid s1 to s2
        nRe = 1;
    else
        error('upscale_geodata:varargin',['Varargin is ' ...
            num2str(varargin{1}) ', which is not a valid option.']);
    end
    
    if numel(varargin(:)) > 1
    	type = varargin{2};
    end
else
    if d1Avg < d2Avg
    	nRe = 1;
    elseif d2Avg < d1Avg %s2 at lower resolution, so upscale to s1
        nRe = 2;
    else %lat and lon grids are not aligned but are at the same resolution
        error('upscale:dSame',['The lat and lon grids are not aligned '...
            'but are at the same resolution.  This case has not been coded for.']);
    end
end
    
if nRe == 1 %s1 at lower resolution, so upscale to s2
    s1 = regrid_geodata(s1,s2,type);
elseif nRe == 2  %s2 at lower resolution, so upscale to s1
    s2 = regrid_geodata(s2,s1,type);
else %lat and lon grids are not aligned but are at the same resolution
    error('upscale:dSame',['The lat and lon grids are not aligned '...
        'but are at the same resolution.  This case has not been coded for.']);
end

end
