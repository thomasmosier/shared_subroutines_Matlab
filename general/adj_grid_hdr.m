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

function hdrNew = adj_grid_hdr(hdrOrg, bndsInd)



hdrNew = nan(6,1);

hdrNew(1) = abs(bndsInd(2) - bndsInd(1)) + 1;
hdrNew(2) = abs(bndsInd(3) - bndsInd(4)) + 1;
hdrNew(3) = hdrOrg(3) + (bndsInd(1) - 1)*hdrOrg(5);
hdrNew(4) = hdrOrg(4) + (hdrOrg(2) - bndsInd(3))*hdrOrg(5);
hdrNew(5) = hdrOrg(5);
hdrNew(6) = hdrOrg(6);

%ESRI ASCII HEADER FORMAT:
%     NCOLS xxx
%     NROWS xxx
%     XLLCENTER xxx | XLLCORNER xxx
%     YLLCENTER xxx | YLLCORNER xxx
%     CELLSIZE xxx
%     NODATA_VALUE xxx