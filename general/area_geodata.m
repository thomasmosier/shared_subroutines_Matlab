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


function area = area_geodata(lon,lat,strTyp)



%Calculate cartesian area of grid over the Earth.

%Inputs: 
%'lat' and 'lon' vectors of latitude and longitude.
%Optional argument indicates whether the input vectors indicate edges ('e') 
%or centroids ('c') of the grid.  Default is centroid.

%Output: 
%'area' = area of each lat/lon cell (units = 'r's units ^2), 

rEarth = 6371000; %units = meters; reference = David R. Lide, ed. Handbook of Chemistry and Physics (81st ed.). CRC. ISBN 0-8493-0481-4.
    

%Ensure that Lat is row vector
if any(size(lat) == 1)
    lat = lat(:);
end
if any(size(lon) == 1)
    lon = lon(:)';
end


if nargin == 2 || isempty(strTyp) || regexpbl(strTyp,'c')
    %Calculate half the difference 
    if numel(lat) > 1
        dLat = 0.5*diff(lat,1,1);
    elseif numel(lat) == 1
        dLat = 90;
    end
    if numel(lon) > 1
        dLon = 0.5*diff(lon,1,2);
    elseif numel(lon) == 1
        dLon = 180;
    end

    %Calculate bounding latitudes:
    if numel(lat) == 1
        latBS = lat - dLat;
        latBN = lat + dLat;
    elseif lat(1) > lat(2)
        latBS = lat - abs([dLat; dLat(end,:)]);
        latBN = lat + abs([dLat(1,:); dLat]);
    elseif lat(1) < lat(2)
        latBS = lat - abs([dLat(1,:); dLat]);
        latBN = lat + abs([dLat; dLat(end,:)]);
    else
       error('area_Geodata:lat','The latitudes in the georeferenced dataset are repeated.') 
    end

%     %Ensure that Lon is column vector
%     [rLon, cLon] = size(lon);
%     if rLon > cLon && (rLon == 1 || cLon == 1)
%         
%     end
    
    %Calculate bounding longitudes:
    if numel(lon) == 1
        lonBW = lon - dLon;
        lonBE = lon + dLon;
    elseif lon(1) < lon(2)
        lonBW = lon - abs([dLon(:,1), dLon]);
        lonBE = lon + abs([dLon, dLon(:,end)]);
    elseif lon(1) > lon(2)
        lonBW = lon - abs([dLon, dLon(:,end)]);
        lonBE = lon + abs([dLon(:,1), dLon]);
    else
       error('area_Geodata:lon','The longitudes in the geo-referenced dataset are repeated.') 
    end
elseif ~isempty(regexpi(strTyp,'e'))
    latBN = lat(1:end-1);
    latBS = lat(2:end);
    lonBW = lon(1:end-1);
    lonBE = lon(2:end);
end

if any(size(lat) == 1) && any(size(lon) == 1) 
    [lonBW,latBS] = meshgrid(lonBW,latBS);
    [lonBE,latBN] = meshgrid(lonBE,latBN);

    %Calculate fractional area over sphere:
    %Specify the reference sphere as Earth with units of 'meters' (output in
    %units of meters^2)
%     area = (2*pi)*rEarth^2*abs(sind(yNE)-sind(ySW)).*(abs(xNE-xSW)/360);
%     area1 = single(areaquad(double(ySW),double(xSW),double(yNE),double(xNE),referenceSphere('earth','meters')));
    
%     area = single(areaquad(double(latBS),double(lonBW),double(latBN),double(lonBE),referenceSphere('earth','meters')));
end

%Built in function from Matlab:
% area = areaquad(latBS,lonBW,latBN,lonBE,referenceSphere('earth','meters'));

%Formula described below:
area = (pi/180)*rEarth^2*abs(sind(latBN)-sind(latBS)).*abs(lonBE-lonBW);

% %Formula based on calculus:
% %***Requires changing latitude reference frame from -90 to +90 to 0 to +180
% area = (pi/180)*rEarth^2*abs(cosd(latBS+90)-cosd(latBN+90)).*abs(lonBE-lonBW);

%DESCRIPTION OF METHOD:
%https://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2004/msg00023.html
% We started with the formula for the area of the earth between a line of 
%latitude and the north pole (the area of a spherical cap, listed in the 
%Dr. Math FAQ on Geometric Formulas).
% A = 2*pi*R*h
% where R is the radius of the earth and h is the perpendicular distance 
%from the plane containing the line of latitude to the pole. We can 
%calculate h using trigonometry as
% h = R*(1-sin(lat))
% Thus the area north of a line of latitude is
% A = 2*pi*R^2(1-sin(lat))
% The area between two lines of latitude is the difference between the area 
%north of one latitude and the area north of the other latitude:
% A = |2*pi*R^2(1-sin(lat2)) - 2*pi*R^2(1-sin(lat1))|
% = 2*pi*R^2 |sin(lat1) - sin(lat2)|
% The area of a lat-lon rectangle is proportional to the difference in the 
%longitudes. The area I just calculated is the area between longitude lines 
%differing by 360 degrees. Therefore the area we seek is
% A = 2*pi*R^2 |sin(lat1)-sin(lat2)| |lon1-lon2|/360
% = (pi/180)R^2 |sin(lat1)-sin(lat2)| |lon1-lon2|

% - Doctor Rick, The Math Forum
% http://mathforum.org/dr.math/ 
