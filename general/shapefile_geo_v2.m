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

function [polyBl, polyVal] = shapefile_geo_v2(S, field, lat, lon)
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
polyVal = polyBl;

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
    keyboard
    [xsplit, ysplit] = polysplit(shpX, shpY);
    isCw = ispolycw(xsplit, ysplit);
    mainPolyIndices = find(isCw);
    nHolesPer = diff([mainPolyIndices;length(isCw)+1]) - 1;

    
    %Loop over points in polygon:
    for ipoly = 1:length(mainPolyIndices)
        isInMain = inpolygon(lonMesh, latMesh, xsplit{mainPolyIndices(ipoly)}, ysplit{mainPolyIndices(ipoly)});
        if nHolesPer(ipoly) > 0
            isInHole = zeros(size(lonMesh), 'single');
            for ihole = 1:nHolesPer(ipoly)
                holeTemp = inpolygon(lonMesh, latMesh, xsplit{mainPolyIndices(ipoly)+ihole}, ysplit{mainPolyIndices(ipoly)+ihole});
                isInHole(holeTemp == 1) = 1;
            end

            polyBl( isInMain & ~isInHole ) = 1;
            polyVal( isInMain & ~isInHole ) = S(zz).(field);
        else
            polyBl(isInMain == 1) = 1; 
            polyVal(isInMain == 1) = S(zz).(field);
        end
    end
end


%
% write_ESRI_v4(polyBl, ESRI_hdr(lon, lat, 'corner'), fullfile(pwd,'UIB_glacier_bl.asc'), 0);

% 
% %Process one million grid points at once
% indEachStep = 10^4;
% 
% %Break up grid into segments:
% %Create indice segmens to loop over:
% indSegments = (1:indEachStep:numel(lonMesh));
% indSegments(end) = numel(lonMesh) + 1;
% 
% %Initialize output:
% if isnumeric(S(1).(attribute)) && isscalar(S(1).(attribute))
%     value = nan(size(lonMesh), 'single');
% else
%     value = cell(size(lonMesh));
% end
% 
% %Enter loop:
% for zz = 1 : numel(indSegments) - 1
%     disp([num2str(zz) ' of ' num2str(numel(indSegments) - 1)]);
%     
%     indCurr = indSegments(zz) : indSegments(zz+1) - 1;
%     
%     SCurr = S;
% %%COMMENTED OUT REMOVING polygons:
% %     %Remove elements from shapefile if significantly outsides bounds:
% %     lonMin = nanmin(nanmin(lonMesh(indCurr))) - bnd;
% %     lonMax = nanmax(nanmax(lonMesh(indCurr))) + bnd;
% %     latMin = nanmin(nanmin(latMesh(indCurr))) - bnd;
% %     latMax = nanmax(nanmax(latMesh(indCurr))) + bnd;
% %     tic
% %     for ii = numel(SCurr) : -1 : 1
% %         if all(SCurr(ii).(strX) < lonMin | isnan(SCurr(ii).(strX))) || all(SCurr(ii).(strX) > lonMax | isnan(SCurr(ii).(strX))) || all(SCurr(ii).(strY) < latMin | isnan(SCurr(ii).(strY))) || all(SCurr(ii).(strY) > latMax | isnan(SCurr(ii).(strY)))
% %             SCurr(ii) = [];
% %         end
% %         if mod(ii,100) == 0
% %            disp([num2str(ii) ' ' num2str(round(toc))]) 
% %         end
% %     end
% 
%     if ~isvector(shpX) || ~isvector(shpY) || length(shpX) ~= length(shpY)
%         error('xv and yv must be vectors of the same length');
%     end
% 
%     %-----------------------------
%     % Find number of starting
%     % indices of polygons
%     %-----------------------------
% 
%     [xsplit, ysplit] = polysplit(shpX, shpY);
%     isCw = ispolycw(xsplit, ysplit);
%     mainPolyIndices = find(isCw);
%     nHolesPer = diff([mainPolyIndices;length(isCw)+1]) - 1;
% 
% %-----------------------------
% % Test if points are in each
% % polygon
% %-----------------------------
% 
% 
%     in = zeros(size(lonMesh), 'single');
%     for ipoly = 1:length(mainPolyIndices)
%         isInMain = inpolygon(lonMesh, latMesh, xsplit{mainPolyIndices(ipoly)}, ysplit{mainPolyIndices(ipoly)});
%         if nHolesPer(ipoly) > 0
%             isInHole = zeros(length(lonMesh), nHolesPer(ipoly), 'single');
%             for ihole = 1:nHolesPer(ipoly)
%                 isInHole(:,ihole) = inpolygon(lonMesh, latMesh, xsplit{mainPolyIndices(ipoly)+ihole}, ysplit{mainPolyIndices(ipoly)+ihole});
%             end
% 
%             crdTemp = isInMain & ~any(isInHole,2);
%         else
%             crdTemp = isInMain;
%         end
%         in(crdTemp == 1) = 1; 
%     end
%     
% 
% %     nclockwise = cellfun(@(x,y) length(find(ispolycw(x,y))), {SCurr.(strX)}, {SCurr.(strY)});
% %     repVals = cellfun(@(nrep, val) repmat({val},1,nrep), num2cell(nclockwise), {SCurr.(attribute)}, 'UniformOutput', false);
% %     repVals = [repVals{:}];
% % 
% %     [~, indexCell] = inpolygons(lonMesh(indCurr), latMesh(indCurr), [SCurr.(strX)], [SCurr.(strY)]);
% % 
% %     try
% %         index = cell2mat(indexCell);
% %     catch
% %         warning('interpshapefile:overlap','Overlapping polygons; using first match for each point');
% %         index = zeros(size(indexCell));
% %         for icell = 1:numel(index)
% %             index(icell) = indexCell{icell}(1);
% %         end
% %     end
% % 
% %     isGood = (index ~= 0);
% %     if isnumeric(SCurr(1).(attribute)) && isscalar(SCurr(1).(attribute))
% %         value(indCurr(isGood)) = cell2mat(repVals(index(isGood)));
% %     else
% %         value(indCurr(isGood)) = repVals(index(isGood));
% %     end
% end