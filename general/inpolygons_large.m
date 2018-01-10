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

function [in, index] = inpolygons_large(x,y,xv,yv)
%INPOLYGONS Performs inpolygon on multiple polygons, holes possible
%
%   in = inpolygons(x,y,xv,yv)
%   [in, index] = inpolygons(x,y,xv,yv)
%
% This function is an extension of Matlab's inpolygon function.  It allows
% the input polygon vertices to describe multiple NaN-delimited polygons.
% The polygons can also include holes.
%
% Input variables:
%
%   x:      x coordinates of points, any dimensions
%
%   y:      y coordinates of points, same dimensions as x
%
%   xv:     x coordinates of polygon vertices, vector.  Separate polygons
%           by NaN.  List vertices of polygons clockwise and holes
%           counterclockwise (use ispolycw to test this).  Holes must
%           immediately follow the polygon with which they are associated.
%
%   yv:     y coordinates of polygon vertices, vector same length and
%           format as xv.  
%
% Output variables:
%
%   in:     matrix in the same size as x and y where in(p,q) = 1 if the
%           point (x(p,q), y(p,q)) is inside any of the polygons defined by
%           xv and yv    
%
%   index:  cell array with the same dimensions as x and y holding the
%           indices of the polygons in which each point was found (0 for
%           point outside all polygons).  
%
% Example:
%
% xv = [1 1 7 7 1 NaN 2 3 3 2 2 NaN 5 6 5 5 NaN 7 8 9 8 7];
% yv = [1 4 4 1 1 NaN 2 2 3 3 2 NaN 2 2 3 2 NaN 8 9 8 7 8];
% x = 10 * rand(20,10); y = 10 * rand(20,10);
% [in, index] = inpolygons(x, y, xv, yv);
% index = cell2mat(index);  % No overlapping polygons allows this.
% [f, v] = poly2fv(xv, yv);
% hold on;
% patch('Faces', f, 'Vertices', v, 'FaceColor', [.9 .9 .9], ...
%       'EdgeColor', 'none');
% plot(x(in), y(in), 'r.', x(~in), y(~in), 'b.');
% plot(x(index==1), y(index==1), 'go', x(index==2), y(index==2), 'mo');


% Copyright 2005 Kelly Kearney

%-----------------------------
% Check inputs
%-----------------------------

if size(x) ~= size(y)
    error('x and y must have the same dimensions');
end

if ~isvector(xv) || ~isvector(yv) || length(xv) ~= length(yv)
    error('xv and yv must be vectors of the same length');
end

%-----------------------------
% Find number of and starting
% indices of polygons
%-----------------------------

[xsplit, ysplit] = polysplit(xv, yv);
isCw = ispolycw(xsplit, ysplit);
mainPolyIndices = find(isCw);
nHolesPer = diff([mainPolyIndices;length(isCw)+1]) - 1;

%-----------------------------
% Test if points are in each
% polygon
%-----------------------------

originalSize = size(x);
x = x(:);
y = y(:);


%UPDATED CODE:
in = zeros(length(x),1);
index = num2cell(zeros(size(x)));
for ipoly = 1:length(mainPolyIndices)
%     if mod(ipoly,100) == 0
%         disp([num2str(ipoly) ' of ' num2str(length(mainPolyIndices)) ' iterations.'])
%     end
    isInMain = inpolygon(x, y, xsplit{mainPolyIndices(ipoly)}, ysplit{mainPolyIndices(ipoly)});
    if nHolesPer(ipoly) > 0
        isInHole = zeros(length(x), nHolesPer(ipoly));
        for ihole = 1:nHolesPer(ipoly)
            isInHole(:,ihole) = inpolygon(x, y, xsplit{mainPolyIndices(ipoly)+ihole}, ysplit{mainPolyIndices(ipoly)+ihole});
        end
        indIn = find(isInMain & ~any(isInHole,2));
    else
        indIn = find(isInMain);
    end
    
    if ~isempty(indIn)
        in(indIn) = 1;
        
        if numel(indIn) == 1
            if index{indIn} == 0
                index{indIn} = ipoly;
            else
                index{indIn} = [index{indIn}, ipoly];
            end
        else
            for iInd = 1 : numel(indIn)
                if index{indIn(iInd)} == 0
                    index{indIn(iInd)} = ipoly;
                else
                    index{indIn(iInd)} = [index{indIn(iInd)}, ipoly];
                end
            end
        end
    end
end

if nargout == 2
    index = reshape(index, originalSize);
end

% %ORIGINAL CODE:
% isIn = zeros(length(x), length(mainPolyIndices));
% for ipoly = 1:length(mainPolyIndices)
%     isInMain = inpolygon(x, y, xsplit{mainPolyIndices(ipoly)}, ysplit{mainPolyIndices(ipoly)});
%     if nHolesPer(ipoly) > 0
%         isInHole = zeros(length(x), nHolesPer(ipoly));
%         for ihole = 1:nHolesPer(ipoly)
%             isInHole(:,ihole) = inpolygon(x, y, xsplit{mainPolyIndices(ipoly)+ihole}, ysplit{mainPolyIndices(ipoly)+ihole});
%         end
%         isIn(:,ipoly) = isInMain & ~any(isInHole,2);
%     else
%         isIn(:,ipoly) = isInMain;
%     end
% end
% 
% in = any(isIn, 2);
% in = reshape(in, originalSize);
% 
% if nargout == 2
% 
%     index = num2cell(zeros(size(x)));
%     for ipoint = 1:length(x)
%         loc = find(isIn(ipoint,:));
%         if ~isempty(loc)
%             index{ipoint} = loc;
%         end
%     end
%     
%     index = reshape(index, originalSize);
% end







