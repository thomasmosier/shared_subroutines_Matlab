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

function dataOut = geodata_area_wgt(lonIn,latIn,dataIn,lonOut,latOut, varargin)



%Varargin specifies 
if ~isempty(varargin(:)) && regexpbl(varargin{1}, 'nansum')
    type = 'nansum';
else
    type = 'sum';
end

%Make sure data does not have superfuluous index:
dataIn = squeeze(dataIn);
if length(size(dataIn)) == 3
   error('area_int:dataD','Function only designed to work with 2D arrays'); 
end

%%Reorient data and reference vectors:
%Transpose reference longitude from row to column vector:
[rLonR, cLonR] = size(lonOut);
if rLonR ~= 1 && cLonR == 1
    lonOut = lonOut';  
elseif rLonR ~= 1 && cLonR ~= 1
    error('area_int_1D:refMat',['The reference longitudes are in matrix '...
        'form but a vector is required.']);
end

%Transpose reference longitude from row to column vector:
[rLonD, cLonD] = size(lonIn);
if rLonD ~= 1 && cLonD == 1
    lonIn = lonIn';  
elseif rLonD ~= 1 && cLonD ~= 1
    error('area_int_1D:refMat',['The reference longitudes are in matrix '...
        'form but a vector is required.']);
end

%Transpose reference latitude from colum to row vector:
[rLatR, cLatR] = size(latOut);
if rLatR == 1 && cLatR ~= 1
    latOut = latOut';  
elseif rLatR ~= 1 && cLatR ~= 1
    error('area_int_1D:refMat',['The reference longitudes are in matrix '...
        'form but a vector is required.']);
end

%Transpose data latitude from colum to row vector:
[rLatD, cLatD] = size(latIn);
if rLatD == 1 && cLatD ~= 1
    latIn = latIn';  
elseif rLatD ~= 1 && cLatD ~= 1
    error('area_int_1D:refMat',['The reference longitudes are in matrix '...
        'form but a vector is required.']);
end

%Flip and sort data so that it is in map orientation:
if ~issorted(lonIn)
    [lonIn, indLonOrdIn] = sort(lonIn);
    dataIn = dataIn(:,indLonOrdIn);
end
if ~issorted(flipud(latIn))
    [latIn, indLatOrdIn] = sort(latIn,'descend');
    dataIn = dataIn(indLatOrdIn,:);
end

%Sort reference coordinates so they're in map orientation:
indLonOrdOut = nan;
if ~issorted(lonOut)
    [lonOut, indLonOrdOut] = sort(lonOut);
end

indLatOrdOut = nan;
if ~issorted(flipud(latOut)) 
    [latOut, indLatOrdOut] = sort(latOut,'descend');
end

%Check that size of lat and lon matches size of data
if length(latIn) ~= length(dataIn(:,1)) || length(lonIn) ~= length(dataIn(1,:))
    error('area_int_2D:dataMismatch',['The dimensions of the data do '...
        'not match the dimensions of the corresponding latitude and '...
        'longitude vectors.']);
end

%Initialize output array:
dataOut = nan(numel(latOut), numel(lonOut));

%If gridding is same, but one is simply outside the other, fill in with
%nans and return
ordrChck = max(-order(nanmean([abs(diff(latOut(:))); abs(diff(lonOut(:)))]))+1,0);

[latTemp, indUseLatOut, indUseLatIn] = intersect(round2(latOut, ordrChck),round2(latIn, ordrChck));
[lonTemp, indUseLonOut, indUseLonIn] = intersect(round2(lonOut, ordrChck),round2(lonIn, ordrChck));
nLatSame = numel(latTemp);
nLonSame = numel(lonTemp);
%Ensure intersected lat indices retain map orientation
if ~issorted(flipud(latOut(indUseLatOut)))
    [~, indTemp] = sort(latOut(indUseLatOut),'descend');
    indUseLatOut = indUseLatOut(indTemp);
end
if ~issorted(flipud(latIn(indUseLatIn)))
    [~, indTemp] = sort(latIn(indUseLatIn),'descend');
    indUseLatIn = indUseLatIn(indTemp);
end

%Create boolean value to indicate grids are actually same but one is contained
%inside the other.  This is much more computationally efficient to solve.
blSame = 0;
if nLatSame == numel(latOut) && nLonSame == numel(lonOut) %Reference grid is subset of input grid
    dataOut = dataIn(indUseLatIn,indUseLonIn);
    blSame = 1;
elseif nLatSame == numel(latIn) && nLonSame == numel(lonIn) %Input grid is subset of reference grid
    dataOut(indUseLatOut,indUseLonOut) = dataIn;
    blSame = 1;
end
    
if blSame == 1
    %Flip output if reference were flipped:
    if ~isnan(indLatOrdOut)
        dataOut = dataOut(indLatOrdOut,:);
    end
    if ~isnan(indLonOrdOut)
        dataOut = dataOut(:,indLonOrdOut);
    end
    
    return
end


%Define grid at edge of cells:
lonOutEdg = box_edg(lonOut);
latOutEdg = box_edg(latOut);
lonInEdg  = box_edg(lonIn);
latInEdg  = box_edg(latIn);

%Approximate formulas being used:
%2*pi*rEarth/360 being dropped because fractional difference is all that
%matters
% dx = 2*pi*rEarth*(dLon/360)*cosd(lat);
% dy = 2*pi*rEarth*(dLat/360);


% %Input distances
% dlLatIn = abs(diff(latInEdg));
% dlLonIn = nan([numel(latIn), numel(lonIn)], 'single');
% dLonIn = abs(diff(lonInEdg));
% for ii = 1 : numel(latIn)
%     %Note: This formulations is a linear approximation
%     dlLonIn(ii,:) = dLonIn*cosd(latIn(ii));
% end
% 
% %Output distances
% dlLatOut = abs(diff(latOutEdg));
% dlLonOut = nan([numel(latOut), numel(lonOut)], 'single');
% dLonOut = abs(diff(lonOutEdg));
% for ii = 1 : numel(latOut)
%     %Note: This formulations is a linear approximation
%     dlLonOut(ii,:) = dLonOut*cosd(latOut(ii));
% end

%Initialize output array:
dataOut = nan([numel(latOut), numel(lonOut)], 'single');

%Find indices to loop over (indices can be skipped if they are entirely outside of input grid cells)
indLonOutUse = find(lonOutEdg(2:end) > min(lonInEdg) & lonOutEdg(1:end-1) < max(lonInEdg));
indLatOutUse = find(latOutEdg(2:end) < max(latInEdg) & latOutEdg(1:end-1) > min(latInEdg));

szDataIn = size(dataIn);

for jj = 1 : numel(indLatOutUse)
    %Bottom of input must be less than top of output & Top of input must be greater than botton of output
    indLatCurr = find(latInEdg(2:end) < latOutEdg(indLatOutUse(jj)) & latInEdg(1:end-1) > latOutEdg(indLatOutUse(jj)+1));
    
    %Find latitudinal overlap:
    latTop = min(latInEdg(indLatCurr  ), latOutEdg(indLatOutUse(jj)  ));
    latBtm = max(latInEdg(indLatCurr+1), latOutEdg(indLatOutUse(jj)+1));
    
    for ii = 1 : numel(indLonOutUse)
        %Right of input must be greater than left of output & Left of input must be less than right of output
        indLonCurr = find(lonInEdg(2:end) > lonOutEdg(indLonOutUse(ii)) & lonInEdg(1:end-1) < lonOutEdg(indLonOutUse(ii)+1));
        
        %Find longitudinal overlap:
        lonLft = max(lonInEdg(indLonCurr  ), lonOutEdg(indLonOutUse(ii)  ));
        lonRgt = min(lonInEdg(indLonCurr+1), lonOutEdg(indLonOutUse(ii)+1));
        
        %Find weighting factors:
        %Area section on sphere in degrees = (4*pi*R^2/360)*(lon2 - lon1)*(cos(lat1)-cos(lat2))
        %Remove constant factor because it factors out in weighting
        %"dly" factor
        dly = (cosd(latBtm) - cosd(latTop))*ones([1,numel(indLonCurr)]);
        %"dlx" factor
        dlx = ones([numel(indLatCurr),1])*(lonRgt - lonLft);
        wgt = dly.*dlx;  
        
        %Create Indices in full form
        rowMesh = indLatCurr*ones([1, numel(indLonCurr)]);
        colMesh = ones([numel(indLatCurr), 1])*indLonCurr;
        
        indCurr = sub2ind(szDataIn, rowMesh(:), colMesh(:));
        
        %Calculate output
        if regexpbl(type, 'nan')
            dataOut(indLatOutUse(jj), indLonOutUse(ii)) = nanwmean(dataIn(indCurr(:)), wgt(:));
        else
            dataOut(indLatOutUse(jj), indLonOutUse(ii)) = wmean(dataIn(indCurr(:)), wgt(:));
        end
    end
end



%Flip output if reference were flipped:
if ~isnan(indLatOrdOut)
    dataOut = dataOut(indLatOrdOut,:);
end
if ~isnan(indLonOrdOut)
    dataOut = dataOut(:,indLonOrdOut);
end

% %Plot output:
% close all
% cSize = 40; %size of circles
% %Pick color bars:
% cMn = min([min(min(data)),min(min(dataR)),min(min(dataC))]);
% cMx = max([max(max(data)),max(max(dataR)),max(max(dataC))]);
% %Plot input data:
% [xD, yD] = meshgrid(lonD,latD);
% figure
% scatter(reshape(xD,1,[]),reshape(yD,1,[]),cSize,reshape(data,1,[]),'fill')
% colorbar
% caxis([cMn, cMx])
% %Plot common array:
% [xC, yC] = meshgrid(lonC,latC);
% figure
% scatter(reshape(xC,1,[]),reshape(yC,1,[]),cSize,reshape(dataC,1,[]),'fill')
% colorbar
% caxis([cMn, cMx])
% % %Plot common array (for reference grid region):
% % [xC, yC] = meshgrid(lonC(indWR:indER-1),latC(indNR:indSR-1));
% % figure
% % scatter(reshape(xC,1,[]),reshape(yC,1,[]),cSize,reshape(dataC(indWR:indER-1,indNR:indSR-1),1,[]),'fill')
% % colorbar
% % caxis([cMn, cMx])
% %Plot reference data:
% [xR, yR] = meshgrid(lonR,latR);
% figure
% scatter(reshape(xR,1,[]),reshape(yR,1,[]),cSize,reshape(dataR,1,[]),'fill')
% %pcolor(xR,yR,dataR)
% colorbar
% caxis([cMn, cMx])
