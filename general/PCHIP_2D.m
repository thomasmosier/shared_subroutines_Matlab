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

function Zo = PCHIP_2D(Xi,Yi,Zi,Xo,Yo)



%Uses Matlab's built in 1D pchip function in 2D.  Input and output formats
%same as for Matlab's 'interp2' function.

warning('OFF','MATLAB:chckxy:IgnoreNaN');

%Ensure input coordinates are correct size
if ~any(size(Xi) == 1) && ~any(size(Yi) == 1) 
    if isequal(size(Xi), size(Yi))
        if Xi(1,1) ==  Xi(2,1)
            Xi = Xi(1,:);
        else
            error('PHCHIP2d:xDiffCrdIn',['The x inputs are a grid, and ' ...
                'values differ between rows. x values cannot differ as a function of the y value).']);
        end
        
        if Yi(1,1) ==  Yi(1,2)
            Yi = Yi(:,1);
        else
            error('PHCHIP2d:yDiffCrdIn',['The y inputs are a grid, and ' ...
                'values differ between rows. y values cannot differ as a function of the x value).']);
        end
    else
        error('PHCHIP2d:xyDiffSizeIn','The x and y inputs are matrices but not of equal shape. This is not allowed.')
    end
end

if ~any(size(Xo) == 1) && ~any(size(Yo) == 1) 
    if isequal(size(Xo), size(Yo))
        if Xo(1,1) ==  Xo(2,1)
            Xo = Xo(1,:);
        else
            error('PHCHIP2d:xDiffCrdOut',['The x output coordinates are a grid, and ' ...
                'values differ between rows. x values cannot differ as a function of the y value).']);
        end
        
        if Yo(1,1) ==  Yo(1,2)
            Yo = Yo(:,1);
        else
            error('PHCHIP2d:yDiffCrdOut',['The y output coordinates are a grid, and ' ...
                'values differ between rows. y values cannot differ as a function of the x value).']);
        end
    else
        error('PHCHIP2d:xyOutDiffSize','The x and y inputs are matrices but not of equal shape. This is not allowed.')
    end
end


%Reform Xa, Ya, Xi, and Yi (Xa to colum vectors and Ya to row vectors):
if length(Xi(:,1)) > 1 && length(Xi(1,:)) == 1 
    Xi = Xi';
end
if length(Xo(:,1)) > 1 && length(Xo(1,:)) == 1 
    Xo = Xo';
end
if length(Yi(1,:)) > 1 && length(Yi(:,1)) == 1 
    Yi = Yi';
end
if length(Yo(1,:)) > 1 && length(Yo(:,1)) == 1 
    Yo = Yo';
end

%Ensure X vectors go from W to E and Y vectors from N to S:
if Xi(1,1) > Xi(1,2)
    Xi = fliplr(Xi);
    Zi = fliplr(Zi);
end
if Yi(1,1) < Yi(2,1)
    Yi = flipud(Yi);
    Zi = flipud(Zi);
end
flipOut = [0,0];
if Xo(1,1) > Xo(1,2)
    Xo = fliplr(Xo);
    flipOut(1) = 1;
end
if Yo(1,1) < Yo(2,1)
    Yo = flipud(Yo);
    flipOut(2) = 1;
end

%Initialize temporary matrix (row ind of original, columns of interp ind):
ZiTemp = nan( length(Xi(:,1)), length(Xo(1,:)) );

for jj = 1 : length( Zi(:,1) )  %Interpolate over each row of input data

    %interpolate, using Matlab's built in 'pchip', which is 1D. 
    %at each of these points and save as row of output data.
    if length(Zi(jj,:)) - sum(isnan(Zi(jj,:))) < 2
        %Cannot perform interpolation with less than two points
        ZiTemp(jj, : ) = NaN; 
        warning('PCHIP_2D:NaN','Interpolation output for current row is NaN because less than two input points in the row are non-NaN.');
    else
        ZiTemp(jj, : ) = pchip( Xi(1,:), Zi(jj,:), Xo(1,:));
    end
end

%Initialize output matrix:
Zo = nan(length(Yo(:,1)), length(Xo(1,:)));

%Interpolate each column of ZiTemp:
for kk = 1 : length( ZiTemp(1,:) )
    if length(ZiTemp(:,kk)) - sum(isnan(ZiTemp(:,kk))) < 2
        Zo( : , kk) = NaN;
        warning('PCHIP_2D:NaN','Interpolation output for current column is NaN because less than two input points in the column are non-NaN.');
    else
        Zo( : , kk) = pchip( Yi(:,1), ZiTemp(:,kk), Yo(:,1) );
    end
end


%%Set points outside grid to nan
Zo = grid_compare_nan(Xi,Yi,Zi,Xo,Yo,Zo);


%If Xi or Yi were changed, flip output:
if flipOut(1) == 1
    Zo = fliplr(Zo);
end
if flipOut(2) == 1
    Zo = flipud(Zo);
end

warning('ON','MATLAB:chckxy:IgnoreNaN');


