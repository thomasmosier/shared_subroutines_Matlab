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

function [cdf, val] = e_cdf(data, varargin)

type = 'sort';
bins = nan;
if ~isempty(varargin(:))
   for ii = 1 : numel(varargin(:))
      if strcmpi(varargin{ii}, 'bins')
          bins = varargin{ii+1};
      elseif strcmpi(varargin{ii}, 'type')
          type = varargin{ii+1};
      end
   end
end
    
if iscell(data)
   error('e_cdf:cell','This function cannot handle cell input.'); 
end

data(isnan(data)) = [];

if ~isnan(bins) && isnumeric(bins) && bins ~= numel(data) && numel(data) > 9
    if regexpbl(type, {'no','sort'}, 'and')
       warning('e_cdf:sort','No sort and bins were both requested, but are not caompatible. Option is being changed to sort.') 
    end
    
    nPtsIn = numel(data);
    
    cdfIn = linspace(0,1,nPtsIn);
    
    cdfTemp = linspace(0,1,nPtsIn*bins);
    
    data = sort(data);
    dataTemp = interp1(cdfIn, data, cdfTemp, 'pchip');
%     dataTemp = interp(double(data(:)),bins);
    data = nan(bins,1);
    cntr = 0;
    for ii = 1 : bins
        data(ii) = nanmean(dataTemp(cntr+1:cntr+nPtsIn));
        cntr = cntr + nPtsIn;
    end
end
    
nData = numel(data);

if numel(data) == 0
   cdf = nan;
   val = nan;
   
   return
end

if strcmpi(type, 'sort')
    cdf = (cumsum(ones(nData,1))) / nData;

    val = sort(data(:),'ascend');
    
    cdf = [0; cdf];
    val = [val(1); val];
elseif regexpbl(type, {'no','sort'}, 'and')
    cdf = nan(numel(data), 1);
    for ii = 1 : numel(data)
        cdf(ii) = sum(data(ii) >= data(:)) / nData;
    end
    val = data(:);
else
    error('eCdf:method',['Method ' type ' is not known.']);
end



