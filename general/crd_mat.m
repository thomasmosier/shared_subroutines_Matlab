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

function [XMat, YMat] = crd_mat(bndsDd, hdr)



%Uses 'bndsDd' which is in the form [lonW, lonE, latS, latN] for the center
%of the outside row/col of data for a geographically located matrix and
%'hdr' which is in ESRI ASCII gridded header format.
%The function produces 'XMat' and 'YMat' corresponding to the output of
%Matlab's 'Meshgrid'.

XVec = bndsDd(1) : hdr(5)  : bndsDd(2);
YVec = bndsDd(4) : -hdr(5) : bndsDd(3);

[XMat, YMat] = meshgrid(XVec, YVec);
end