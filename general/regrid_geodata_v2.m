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

function sData = regrid_geodata_v2(sData,sRef, varUse, varargin)


%Regrids data in sData to the grid in sRef using an area-integration scheme
%Inputs:
%'sData' is structure array (Downscaling Package format) to be regridded
%'sRef' is structure array (Downscaling Package format) which 'SData' is
%regridded to.
%Outputs:
%'sData.(varOut)' are data of sData regridded to sRef
%'sData.(varLat)Re' are latitudes of regridded reference structure
%'sData.(varLon)Re' are longitudes of regridded reference structure


varLon = 'longitude';
varLat = 'latitude';

if ~isempty(varargin(:))
   type = varargin{1};
else
    type = 'nansum';
end
    
    
varOut = [varUse 'Re'];

%Regrid data:
if length(size(sData.(varUse))) == 3 %If input has time dimension:
    sData.(varOut) = nan(length(sData.(varUse)(:,1,1)),length(sRef.(varLat)),length(sRef.(varLon)),'single');
    
    for ii = 1 : length(squeeze(sData.(varUse)(:,1,1)))
        sData.(varOut)(ii,:,:) = area_int_2D_v2(sData.(varLon),sData.(varLat),squeeze(sData.(varUse)(ii,:,:)),sRef.(varLon),sRef.(varLat), type);
    end
elseif length(size(sData.(varUse))) == 2 %If input data has no time dimension:
    sData.(varOut) = area_int_2D_v2(sData.(varLon),sData.(varLat),sData.(varUse),sRef.(varLon),sRef.(varLat), type);
end

%Add field to SData with lats and lons
sData.([varLon 'Re']) = sRef.(varLon);
sData.([varLat 'Re']) = sRef.(varLat);
