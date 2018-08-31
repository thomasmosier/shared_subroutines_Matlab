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

function [polyBl, varargout] = shapefile_geo_v2(S, lat, lon, blHole, varargin)
%Note: same as 'interpshapefile' except that it removes shapes that are
%sufficiently outside of lat/lon inputs before doing interpolation. Also
%breaks interpolation up into pieces to reduce memory cost.

%INTERPSHAPEFILE Interpolates data based on a shapefile data structure
%
% value = interpshapefile(S, lat, lon, attribute)
%
% This program determines the attribute value at the given points, based on
% the data in the geographic data structure S (created via
% shaperead).  
%
% If the specified attribute holds a numeric scalar, value will be a
% numeric array; points located outside the polygons will receive a NaN.
% Otherwise, value will be a cell array; points outside the polygons will
% hold an empty array.
%
% NOTE: This program assumes that the geographic region covered by two
% adjacent vertices is small enough that the lines connecting the polygon
% vertices can be approximated using straight lines in a cartesian plane,
% and that the region does not cross the poles.
%
% Input variables:
%   
%   S:              geographic data structure returned by a shaperead.  At
%                   minimum, the structure must have the Lat, Lon, and
%                   specified attribute fields.  If the Geometry field is
%                   present, the function uses only those elements
%                   designated Polygons; if the field is not present (e.g.
%                   if it was deleted in previous processing), it assumes
%                   all elements are polygons.
%
%   lat:            latitude coordinates of points to interpolate
%
%   lon:            longitude coordinates of points to interpolate
%
%   attribute:      name of attribute field to be interpolated
%
% Output variables:
%
%   value:          array of interpolated values corresponding to each
%                   lat/lon input.  If the specified attribute holds a
%                   numeric scalar, value will be a numeric array; points
%                   located outside the polygons will receive a NaN.
%                   Otherwise, value will be a cell array; points outside
%                   the polygons will hold an empty array.    
%
% Example:
% 
%   States = shaperead('usastatelo', 'UseGeoCoords', true);
%   latcorners = [max([States.Lat]) min([States.Lat])];
%   loncorners = [min([States.Lon]) max([States.Lon])];
%   [lat, lon] = track(latcorners, loncorners);
%   statename = interpshapefile(States, lat, lon, 'Name');
%   
%   axesm('mercator');
%   geoshow('usastatelo.shp', 'FaceColor', [.5 .5 1]);
%   plotm(lat, lon, 'r.');
%   textm(lat, lon, statename);
%
% See also:
%
%   shaperead

% Copyright 2005 Kelly Kearney

%------------------------------
% Check input
%------------------------------

% if size(lat) ~= size(lon)
%     error('x and y must have the same dimensions');
% end

%Check if shape structure has 'lon'/'lat' or 'X'/'Y' fields
if all(isfield(S, {'lon','lat'}))
    strX = 'lon';
    strY = 'lat';
elseif all(isfield(S, {'X','Y'}))
    strX = 'X';
    strY = 'Y';
elseif all(isfield(S, {'x','y'}))
    strX = 'x';
    strY = 'y';
else
    error('interpshapefile:unknownCrd','The coordinate fields could not be found.');
end


if ~all(size(lat) == size(lon)) || (any(size(lat) == 1) && any(size(lon) == 1))
    [lonMesh, latMesh] = meshgrid(lon,lat);
elseif ~all(size(lat) == size(lon))
    error('interpshapefile:latLon','lat and lon must be matrices of the same size');
else
    lonMesh = lon;
    latMesh = lat;
end


%%COMMENTED OUT REMOVING polygons:
% %Remove elements from shapefile if significantly outsides bounds:
% lonMin = nanmin(nanmin(lonMesh)) - bnd;
% lonMax = nanmax(nanmax(lonMesh)) + bnd;
% latMin = nanmin(nanmin(latMesh)) - bnd;
% latMax = nanmax(nanmax(latMesh)) + bnd;
% for ii = numel(S) : -1 : 1
%     if all(S(ii).(strX) < lonMin | isnan(S(ii).(strX))) || all(S(ii).(strX) > lonMax | isnan(S(ii).(strX))) || all(S(ii).(strY) < latMin | isnan(S(ii).(strY))) || all(S(ii).(strY) > latMax | isnan(S(ii).(strY)))
%         S(ii) = [];
%     end
% end
    
%------------------------------
% Remove any non-polygon 
% objects
%------------------------------


% If Geometry field is not present, all elements are assumed to be
% polygons.

if isfield(S, 'Geometry')
    ispolygon = strcmp('Polygon', {S.Geometry});
    S(~ispolygon) = [];
end

%------------------------------
% Find attribute value at each
% lat/lon point
%------------------------------


%-----------------------------
% Test if points are in each
% polygon
%-----------------------------

%Initialize output:
polyBl = zeros(size(lonMesh), 'single');
if ~isempty(varargin)
    flds = varargin;
    polyVal = cell(numel(flds(:)), 1);
    [polyVal{:}] = deal(nan(size(lonMesh), 'single'));
else
    polyVal = [];
    if nargout > 1
       error('shapefileGeo:NoInputFields','Field outputs requested, but no input field names provided.') 
    end
end

%Test whether main polygons are CW or CCW orientation:
%In general: CW (1) = main polygon; CCW (0) = hole
%Sometimes the reverse is true

nLp = min([length(S), 50]);
testCw = nan(nLp, 1);

for zz = 1 : length(S) 
    shpX = S(zz).(strX);
    shpY = S(zz).(strY);

    if ~isvector(shpX) || ~isvector(shpY) || length(shpX) ~= length(shpY)
        error('xv and yv must be vectors of the same length');
    end

    %-----------------------------
    % Find number of starting
    % indices of polygons
    %-----------------------------
    [xsplit, ysplit] = poly_split(shpX, shpY);
	isCw = ispoly_cw(xsplit, ysplit);
    
    testCw(zz) = round(nanmean(isCw));
end
%Empirically determine if typical orientation is CW or CCW (use in loop below):
blMainCW = round(nanmean(testCw));


%Loop over polygons:
for zz = 1 : length(S)
%     if mod(zz, 200) == 0
%        disp([num2str(zz) ' of ' num2str(length(S)) ' ' num2str(round(toc)) ]); 
%     end
        
    shpX = S(zz).(strX);
    shpY = S(zz).(strY);

    if ~isvector(shpX) || ~isvector(shpY) || length(shpX) ~= length(shpY)
        error('xv and yv must be vectors of the same length');
    end

    %-----------------------------
    % Find number of starting
    % indices of polygons
    %-----------------------------
    [xsplit, ysplit] = poly_split(shpX, shpY);
	isCw = ispoly_cw(xsplit, ysplit); 
    

    if blMainCW
        mainPoly = find(isCw);
    else
        mainPoly = find(~isCw);
    end
    holePoly = setdiff((1:numel(isCw)), mainPoly);
    nMain = numel(mainPoly);
    nHole = numel(holePoly);
    
    %Loop over features in polygon:
    isInMain = zeros(size(lonMesh), 'single');
    if nMain > 0
        for xx = 1 : nMain
            mainTemp = inpolygon(lonMesh, latMesh, xsplit{mainPoly(xx)}, ysplit{mainPoly(xx)});
            isInMain(mainTemp == 1) = 1;
        end
    end
        
    isInHole = zeros(size(lonMesh), 'single');
    if nHole > 0 && blHole
        for xx = 1 : nHole
            holeTemp = inpolygon(lonMesh, latMesh, xsplit{holePoly(xx)}, ysplit{holePoly(xx)});
            isInHole(holeTemp == 1) = 1;
        end
    end

    %Assign boolean values
    polyBl(isInMain & ~isInHole) = 1;
    
    %Assign non-boolean values based on shapefile fields
    if ~isempty(polyVal)
        for hh = 1 : numel(flds(:))
            if isfield(S(zz), flds{hh}) && ~isempty(S(zz).(flds{hh}))
                valCurr = S(zz).(flds{hh});
            else
                valCurr = 1;
            end
            polyVal{hh}(isInMain & ~isInHole) = valCurr;
        end
    end
end

%Variable output argument:
if nargout > 1
    if nargout - 1 == numel(flds(:))
        varargout = polyVal;
    else
        error('shapefileGeo:inputOutputMismatch', ['There are ' num2str(nargout-1) ' output arguments for fields but ' num2str(numel(flds(:))) ' fields.']);
    end
end

