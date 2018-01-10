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

function dataOut = area_int_2D_v2(lonIn,latIn,dataIn,lonOut,latOut, varargin)

warning('area_int_2D:obsolete','This function is obsolete. Use geodata_area_wgt() instead.')
%Fields used in input data:
%data
%latD
%lonD
%areaD
%latR
%lonR
%areaR
%varargin{1} = 'nansum' (nan's do not affect processing), 'sum' (output is
%nan if any input is nan

%Output:
%'dataOut' = input data that has been area-weighted to output grid 


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
    [lonIn, indLonD] = sort(lonIn);
    dataIn = dataIn(:,indLonD);
end
if ~issorted(flipud(latIn))
    [latIn, indLatD] = sort(latIn,'descend');
    dataIn = dataIn(indLatD,:);
end

%Sort reference coordinates so they're in map orientation:
indLonR = nan;
if ~issorted(lonOut)
    [lonOut, indLonR] = sort(lonOut);
end

indLatR = nan;
if ~issorted(flipud(latOut)) 
    [latOut, indLatR] = sort(latOut,'descend');
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
    if ~isnan(indLatR)
        dataOut = dataOut(indLatR,:);
    end
    if ~isnan(indLonR)
        dataOut = dataOut(:,indLonR);
    end
    
    return
end


%Define grid at edge of cells:
lonOutEdg = box_edg(lonOut);
latOutEdg = box_edg(latOut);
lonInEdg  = box_edg(lonIn);
latInEdg  = box_edg(latIn);


%%Create common data grid:
latCEdg = sort(unique([latInEdg;latOutEdg]),'descend');
dLatC = abs(diff(latCEdg));
latC = latCEdg(1:end-1)-0.5*dLatC;

lonCEdg = sort(unique([lonInEdg,lonOutEdg]));
dLonC = abs(diff(lonCEdg));
lonC = lonCEdg(1:end-1)-0.5*dLonC;

%Initialize common data array and output array:
dataC = nan(length(latC),length(lonC));

%%Integrate data from common grid to reference grid:
%Calculate areas of grids:
areaIn = area_geodata(lonInEdg,latInEdg,'e');
areaIn(isnan(dataIn)) = nan;
areaOut = area_geodata(lonOutEdg,latOutEdg,'e');

areaC = area_geodata(lonCEdg,latCEdg,'e');
%%Place all data within common grid:
for ii = 1 : length(latC)
    for jj = 1 : length(lonC)
        indLatC = find(latInEdg(1:end-1) >= latC(ii) & latInEdg(2:end) <= latC(ii) );

        indLonC = find(lonInEdg(2:end) >= lonC(jj) & lonInEdg(1:end-1) <= lonC(jj) );
        
        if ~isempty(indLatC) && ~isempty(indLonC)
            if length(indLatC) == 2 || length(indLonC) == 2
                dataC(ii,jj) = sum2d(areaIn(indLatC, indLonC).*dataIn(indLatC, indLonC))./sum2d(areaIn(indLatC, indLonC));
            else
                dataC(ii,jj) = dataIn(indLatC, indLonC);
            end
        end
    end
end




areaC(isnan(dataC)) = nan;


%Find all data within each reference grid cell:
for ii = 1 : length(latOut)
    for jj = 1 : length(lonOut)
        indWR = find( round2(lonCEdg,3) == round2(lonOutEdg(jj)  ,3));
        indER = find( round2(lonCEdg,3) == round2(lonOutEdg(jj+1),3));
        indNR = find( round2(latCEdg,3) == round2(latOutEdg(ii)  ,3));
        indSR = find( round2(latCEdg,3) == round2(latOutEdg(ii+1),3));

        if ~isempty(indWR) && ~isempty(indER) && ~isempty(indNR) && ~isempty(indSR)
            if regexpbl(type,'nansum')
                dataOut(ii,jj) = nansum(nansum(areaC(indNR:indSR-1,indWR:indER-1).*dataC(indNR:indSR-1,indWR:indER-1) )) ./ nansum(nansum(areaC(indNR:indSR-1,indWR:indER-1)));
            else %If one or more cells of dataC are NaN, output will be NaN
                dataOut(ii,jj) = sum(sum(areaC(indNR:indSR-1,indWR:indER-1).*dataC(indNR:indSR-1,indWR:indER-1) )) ./ sum(sum(areaC(indNR:indSR-1,indWR:indER-1)));
            end
%             if regexpbl(type,'nansum')
%                 dataR(ii,jj) = nansum(nansum(areaC(indNR:indSR-1,indWR:indER-1).*dataC(indNR:indSR-1,indWR:indER-1) )) ./ areaR(ii,jj);
%             else %If one or more cells of dataC are NaN, output will be NaN
%                 dataR(ii,jj) = sum(sum(areaC(indNR:indSR-1,indWR:indER-1).*dataC(indNR:indSR-1,indWR:indER-1) )) ./ areaR(ii,jj);
%             end
        else
            error('area_int_2D:noPtsC2R',['For ii = ' num2str(ii) ' and jj = ' num2str(jj) ', no points were selected.'])
        end
    end
end

%If areaR is 0, set dataR to 0 (may occur if latitude extends beyond +90 or
%-90)
dataOut(areaOut == 0) = 0;


%Flip output if reference were flipped:
if ~isnan(indLatR)
    dataOut = dataOut(indLatR,:);
end
if ~isnan(indLonR)
    dataOut = dataOut(:,indLonR);
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
