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

function [n, c] = ecdf_hist( data , nBin)
%n = number of elements in each bin
%c = nanmean of all values in each bin

if iscell(data)
   error('e_cdf:cell','This function cannot handle cell input.'); 
end

[cdf, val] = e_cdf( data );
cdf = cdf(2:end);
val = val(2:end);

cdfCent = linspace(0,1,nBin+1);

c = nan(nBin,1);
n = nan(nBin,1);
for ii = 1 : nBin
    ind = find(cdf > cdfCent(ii) & cdf < cdfCent(ii+1));
    c(ii) = nanmean(val(c));
    n(ii) = numel(ind);
end