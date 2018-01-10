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

function [pJoint, xC, yC] = joint_hist(x,y,varargin)


if numel(varargin) == 1
    nBin = varargin{1};
    
    [~, valX] = e_cdf(x, 'bins', nBin, 'type', 'no_sort');
    [~, valY] = e_cdf(y, 'bins', nBin, 'type', 'no_sort');
elseif numel(varargin) == 2
    valX =  varargin{1};
    valY =  varargin{2};
    nBin = numel(valX);
else
    error('joint_hist:arguments','The number of inputs is not acceptable.')
end



if all(isnan(valX)) || all(isnan(valY)) 
    pJoint = nan;
    xC = nan;
    yC = nan;
    return
elseif valX(1) == valX(2) && valY(1) == valY(2)
    valX = valX(2:end);
    valY = valY(2:end);
    nBin = nBin - 1;
end

pJoint = zeros(nBin);
xC = nan(nBin);
yC = nan(nBin);

for ii = 1 : nBin %Loop over x
    if ii == 1
        indX = find(x <= valX(ii));
    else
        indX = find(x > valX(ii-1) & x <= valX(ii));
    end
    
    for jj = 1 : nBin %Loop over y
        if jj == 1
            indY = find(y(indX) <= valY(jj));
        else
            indY = find(y(indX) > valY(jj-1) & y(indX) <= valY(jj));
        end
        
        indJ = indX(indY);
        if ~isempty(indJ)
            pJoint(ii,jj) = numel(indJ);
            xC(ii,jj) = nanmean(x(indJ));
            yC(ii,jj) = nanmean(y(indJ));
        end
    end
end

%Normalize nJoint:
pJoint = 100*pJoint/sum2d(pJoint);

