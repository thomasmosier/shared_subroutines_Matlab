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

function [dataTs, hdr, metaESRI, dataNorms, dataNStns] = read_GPCC_ASCII(fileName)



%Reads GPCC ASCII data in the following format:
    % 14       : rows to skip
    % GPCC Full Data Product Version 006 , produced december 2011
    % (GAUGE ANOMALY-ANALYSIS SPHEREMAP + GPCC CLIMAT v2011)
    % =====================================================
    % Info     : gpcc.dwd.de --  gpcc@dwd.de
    % =====================================================
    % Grid     : 0.50 degree
    % Area     : LON -180.00 to 180.00 LAT 90.00 to -90.00
    % Month    : 01 1950
    % =====================================================
    % column 1 : precipitation totals in mm/month
    % column 2 : gpcc normals v2011 in mm/month
    % column 3 : number of gauges per grid
    % =====================================================
%where the rows begin at 90deg N, 180deg W and proceed eastward.  After
%reaching 90deg N, 180deg E, the next element is for the row immediately
%south of the first, etc.

%The output hdr is in ESRI Hdr format: 
    %line 1 = number of columns
    %line 2 = number of rows
    %line 3 = longitude, lower left corner (center)
    %line 4 = latitude, lower left corner (center)
    %line 5 = cellsize
    %line 6 = no data value

if isunix
    fileName = [filesep fileName];
end

fidData = fopen(fileName, 'r');
     
[hdrLn, ~] = fscanf(fidData, '%d %*s \n',1);

% Skip header lines
for ii = 1: hdrLn
    fgetl(fidData);
end

[tempData, ~] = fscanf(fidData, '%f %f %f \n',[3,inf]);
    
tempData = tempData';

fclose(fidData);

%Reshape data:
dataTs    = reshape(tempData(:,1), 720, 360)';
dataNorms = reshape(tempData(:,2), 720, 360)';
dataNStns = reshape(tempData(:,3), 720, 360)';

dataTs(dataTs == -99999.99) = NaN;
dataNorms(dataNorms == -99999.99) = NaN;
dataNStns(dataNStns == -99999.99) = NaN;

hdr = [720; 360; -180; -90; 0.5; -9999];
metaESRI = {'NCOLS';'NROWS';'XLLCORNER';'YLLCORNER';'CELLSIZE';'NODATA_VALUE'};

dataTs(      dataTs == hdr(6)) = NaN;
dataNorms(dataNorms == hdr(6)) = NaN;
dataNStns(dataNStns == hdr(6)) = NaN;

dataTs = single(dataTs);
end