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

function [lat, lon] = ESRI_hdr2geo(hdr,meta)


%'lat' and 'lon' are centroids of grid boxes defined by inputs 'hdr' and
%'meta', which are 6 row ESRI formatted headers.

if sum2d(~isnan(hdr)) == 0
    lat = nan;
    lon = nan;
    return
end

if ~isempty(regexpi(meta{4},'yllcorner'))
    lat = (hdr(4) + (hdr(2)-0.5)*hdr(5) : -hdr(5) : hdr(4) + 0.5*hdr(5))';
elseif ~isempty(regexpi(meta{4},'yllcenter'))
    lat = (hdr(4) + (hdr(2)-1)*hdr(5) : -hdr(5) : hdr(4))';
else
    error('ESRI_hdr2geo:UnkownRefLat',['The latitude reference is ' meta{4} ', which is an unkown option.']);
end
sData.attLat = {'long_name', 'latitude'; 'units', 'degrees_north'};

%Write lon:
if ~isempty(regexpi(meta{3},'xllcorner'))
    lon = (hdr(3) + 0.5*hdr(5) : hdr(5) : hdr(3) + (hdr(1)-0.5)*hdr(5));
elseif ~isempty(regexpi(meta{3},'xllcenter'))
    lon = (hdr(3) : hdr(5) : hdr(3) + (hdr(1)-1)*hdr(5));
else
    error('ESRI_hdr2geo:UnkownRefLon',['The longitude reference is ' meta{3} ', which is an unkown option.']);
end

