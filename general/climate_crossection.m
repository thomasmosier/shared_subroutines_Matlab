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

function [transectData, hdr] = climate_crossection(climateStr, transectDd)
crop = 'out';   %Out is preferable for longitude, but maybe not for latitude

if regexpi(climateStr,'.asc')
    [transectData, hdr] = read_ESRI(climateStr);
elseif regexpi(climateStr,'.bil')
    %Read and load header of WorldClim binary file:
    [pathClimate, nameClimate, extClimate] = fileparts(climateStr); 
    hdr = read_wc_bin_hdr(fullfile(pathClimate, [nameClimate '.hdr']));    
end

transectDd = [transectDd(1), transectDd(2), (transectDd(3)-hdr(5)),(transectDd(3)+hdr(5))];

if regexpi(climateStr,'.bil')
    %Load WC binary file:
    [transectData, hdr, ~, ~] = read_wc_bin_data(climateStr, hdr, transectDd, crop);  
elseif regexpi(climateStr,'.asc')
    [bndsDd, bndsInd] = adjust_bounds(hdr, transectDd, crop);
    
    if bndsInd(1) > hdr(1) || bndsInd(2) > hdr(1)
        error('The transect longitudinal bounds extend beyond the geographic scope of the file chosen.');
    elseif bndsInd(3) > hdr(2) || bndsInd(4) > hdr(2)
        error(['The transect' char(39) 's latitude is not within the bounds of one of the files chosen.']);
    end

    transectData = transectData(bndsInd(4) : bndsInd(3) , bndsInd(1) : bndsInd(2));
    transectData(transectData == hdr(6)) = NaN;
    
    hdr(1) = length(transectData(1,:));
    hdr(2) = length(transectData(:,1));
    hdr(3) = bndsDd(1) - hdr(5)/2;
    hdr(4) = bndsDd(3) - hdr(5)/2;
end


end
