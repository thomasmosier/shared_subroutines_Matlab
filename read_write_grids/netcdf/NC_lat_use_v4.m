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

function [gcmLat, gcmLatUseInd] = NC_lat_use_v4(gcmLat, attLat, sMeta)


%Find indices for latitude of interest in NetCDF or GRB file:
latUnitInd = strcmpi(attLat,'units');
[latRow, latCol] = find(latUnitInd == 1);

%Find lat edges of GCM grid:
latEdg = box_edg(gcmLat);

%Check units:
latInfo = attLat{latRow, latCol+1};

%Determine frame to use (around coordinates bounds)
if isfield(sMeta,'frame') 
    frame = sMeta.frame;
else
    frame = 0; 
end

if isfield(sMeta,'crop') 
    cropType = sMeta.crop;
elseif isfield(sMeta,'cropType') 
    cropType = sMeta.cropType;
else
    cropType = 'out'; 
end


if isfield(sMeta,'crd') && all(~isnan(sMeta.crd(3:4)))
    latMin = min(sMeta.crd(3:4));
    latMax = max(sMeta.crd(3:4));
    
    %Find indices of file corresponding to desired region:
    if ~isempty(regexpi(latInfo,'deg')) && ~isempty(regexpi(latInfo,'north'))
        if gcmLat(2) > gcmLat(1) %Top of array is South
            
            warning('latUse:cropIndices', 'Check cropping indices for this case');
            if strcmpi(cropType, 'in')
                indN = find(latEdg(2:end) <= latMax, 1, 'last');
                indS = find(latEdg(1:end-1) >= latMin, 1, 'first');
            elseif strcmpi(cropType, 'out')
                indN = find(latEdg(2:end) >= latMax, 1, 'first');
                indS = find(latEdg(1:end-1) <= latMin, 1, 'last');
            else
                error('cro_geodata:unknowntype', ['Method type ' cropType ' has not been programmed for.']);
            end
            
%             latEdg(indS)
%             latEdg(indN+1)
            
            if isempty(indS)
                indS = 1; 
            end
            if isempty(indN)
                indN = numel(gcmLat); 
            end
            
            indN = indN + frame;
            indS = indS - frame;
            
            if indS < 1
                indS = 1;
                warning('crop_geostrucct:outLatN','The Northern latitude cropping indice is out of bounds and therefore being set to 1.');
            end
            if indN > numel(gcmLat)
                indN = numel(gcmLat);
                warning('crop_geostrucct:outLatN','The Northern latitude cropping indice is out of bounds and therefore being set to 1.');
            end
            
            gcmLatUseInd = (indS : indN);
            
            
            
%             if strcmpi(cropType, 'in')
%                 indN = find(latEdg <= latMax, 1, 'first');
%                 indS = find(latMin <= latEdg, 1, 'last');
%             elseif strcmpi(cropType, 'out')
%                 indN = find(latEdg >= latMax, 1, 'last');
%                 indS = find(latMin >= latEdg, 1, 'first');
%             else
%                 error('cro_geodata:unknowntype', ['Method type ' cropType ' has not been programmed for.']);
%             end
%             
%             if isempty(indS)
%                 indS = 1; 
%             end
%             if isempty(indN)
%                 indN = numel(gcmLat); 
%             end
%             
%             %Find latitude indices that are within the desired region:
%             indS = find(latEdg(1:end-1) <= sMeta.crd(3), 1, 'last');
%             if isempty(indS)
%                 indS = 1; 
%             end
%             indN = find(latEdg(2:end) >= sMeta.crd(4), 1, 'first');
%             if isempty(indN)
%                 indN = numel(gcmLat); 
%             end

%             %Indices of GCM cells to use, including double wide border
%             gcmLatUseInd = (indS - frame : indN - 1 + frame); %Empirically this provides two centroids outside of requested bounds
        elseif gcmLat(1) > gcmLat(2) %Top of array is North
            if strcmpi(cropType, 'in')
                indN = find(latEdg(1:end-1) <= latMax, 1, 'first');
                indS = find(latMin <= latEdg(2:end), 1, 'last');
            elseif strcmpi(cropType, 'out')
                indN = find(latEdg(1:end-1) >= latMax, 1, 'last');
                indS = find(latMin >= latEdg(2:end), 1, 'first');
            else
                error('cro_geodata:unknowntype', ['Method type ' cropType ' has not been programmed for.']);
            end
            
            if isempty(indN)
                indN = 1; 
            end
            if isempty(indS)
                indS = numel(gcmLat); 
            end
            
            indN = indN - frame;
            indS = indS + frame;
            
            if indN < 1
                indN = 1;
                warning('crop_geostrucct:outLatN','The Northern latitude cropping indice is out of bounds and therefore being set to 1.');
            end
            if indS > numel(gcmLat)
                indS = numel(gcmLat);
                warning('crop_geostrucct:outLatN','The Northern latitude cropping indice is out of bounds and therefore being set to 1.');
            end
            
            gcmLatUseInd = (indN : indS);
%             %Find latitude indices that are within the desired region:
%             indN = find(latEdg(1:end-1) >= sMeta.crd(4), 1, 'last');
%             if isempty(indN)
%                 indN = 1; 
%             end
%             indS = find(latEdg(2:end) <= sMeta.crd(3), 1, 'first');
%             if isempty(indS)
%                 indS = numel(gcmLat); 
%             end
            
%             %Indices of GCM cells to use, including border
%             if frame == 0
%                 gcmLatUseInd = (indN : indS);
%             else
%                 gcmLatUseInd = (indN + 1 - frame : indS + frame); %Empirically this provides frame centroids outside of requested bounds
%             end
            
        else
            error('NC_lat_use:noChange',['The spacing between latitudes '...
                'in the current NetCDF dataset appears to be 0.  This case has '...
                'not been programmed.']);
        end

    else %No other conditionals written for different units
        error('NC_lat_use:coordSys',['The latitudinal coordinate system '...
            'in the current NetCDF dataset is ' latUnits ', which has not '...
            'been coded for.']);
    end
else
    gcmLatUseInd = (1:numel(gcmLat));
end

%Remove any latitudes above or below North or South Poles, respectively.
indLes = find(gcmLatUseInd < 1);
if ~isempty(indLes)
    gcmLatUseInd(indLes) = 1;
end
indGre = find(gcmLatUseInd > numel(gcmLat));
if ~isempty(indGre)
    gcmLatUseInd(indGre) = numel(gcmLat);
end
%Keep only unique values:
gcmLatUseInd = unique(gcmLatUseInd,'stable');
%List latitudes within target region:
gcmLat = gcmLat(gcmLatUseInd); 
