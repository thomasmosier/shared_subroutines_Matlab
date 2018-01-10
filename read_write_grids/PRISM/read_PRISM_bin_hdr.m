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

function [hdrEsri, metaEsri, bits] = read_PRISM_bin_hdr(pathHdr)

fidHdr = fopen(pathHdr,'r'); 
hdrRaw = textscan(fidHdr,'%s');
fclose(fidHdr);

%Information in PRISM header:
%     'BYTEORDER'
%     'I'
%     'LAYOUT'
%     'BIL'
%     'NROWS'
%     '3105'
%     'NCOLS'
%     '7025'
%     'NBANDS'
%     '1'
%     'NBITS'
%     '32'
%     'BANDROWBYTES'
%     '28100'
%     'TOTALROWBYTES'
%     '28100'
%     'PIXELTYPE'
%     'FLOAT'
%     'ULXMAP'
%     '-125.016666666667'
%     'ULYMAP'
%     '49.9333333333323'
%     'XDIM'
%     '0.008333333333333'
%     'YDIM'
%     '0.008333333333333'
%     'NODATA'
%     '-9999'

bits = str2double(hdrRaw{1}(12));

ncol = str2double(hdrRaw{1}(8));
nrow = str2double(hdrRaw{1}(6));

ulx = str2double(hdrRaw{1}(20));
uly = str2double(hdrRaw{1}(22));

xStep = str2double(hdrRaw{1}(24));
yStep = str2double(hdrRaw{1}(26));

noData = str2double(hdrRaw{1}(28));

if xStep ~= yStep
   error('read_PRISM_bin_hdr:steps', 'The x and y step sizes in the PRISM data are different. This has not been programmed for.') 
end

llx = ulx - 0.5*xStep;
lly = uly - (nrow-0.5)*yStep;

%Write ESRI format header:
%NCOLS      xxx
%NROWS      xxx
%XLLCORNER  xxx
%YLLCORNER  xxx
%CELLSIZE   xxx
%NODATA_VALUE

metaEsri = {'NCOLS';'NROWS';'XLLCORNER';'YLLCORNER';'CELLSIZE';'NODATA_VALUE'};

hdrEsri = [ ncol; ...     %ncols
            nrow; ...     %nrows
            llx; ...  %xllcorner
            lly; ...    %yllcorner
            xStep; ...  %cellsize
            noData];    %NODATA_value
