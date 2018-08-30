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

function [prob, valMid, varargout] = e_pdf(data, nBins, varargin)
%An empirical PDF is essentially a histogram

%Variable input arguments:
rng = nan;
typ = 'frac';
noVal = 'zero';
if ~isempty(varargin)
   for ii = 1 : numel(varargin(:))
      if regexpbl(varargin{ii}, {'rng', 'range'})
          rng = varargin{ii+1};
          if numel(rng) ~= 2
              error('ePdf:rngNotValidSize', 'Range variable input argument is not of length 2.')
          end
      elseif regexpbl(varargin{ii}, 'type')
        typ = varargin{ii+1};
      elseif regexpbl(varargin{ii}, {'no', 'value'})
        noVal = varargin{ii+1};
      end
   end
end


if numel(rng) <= 2 && all(isnan(rng))
    rng = [nanmin(data), nanmax(data)];
end


%Ensure nBins is valid and integer
if isnan(nBins) || ~isnumeric(nBins) || (isnumeric(nBins) && nBins < 1)
   error('ePdf:bins','An unrecognizable bin input wsa provided. It must be a positive integer.'); 
end
nBins = round(nBins);

nPtUse = numel(data);
% nPtUse = sum(data >= rng(1) & data <= rng(2));
% if nPtUse == 0
%     pdf = nan(nBins,1);
%     val = pdf;
%     return
% end

%Create bins:
bins = linspace(rng(1), rng(2), nBins+1);

%Calculate PDF
if regexpbl(noVal, 'zero')
    prob = zeros(nBins,1);
elseif regexpbl(noVal, 'nan')
    prob = nan(nBins,1);
else
    error('ePdf:noVal', ['The option for what to do when the data has no value is ' noVal ', which is not an allowed option.']);
end
valAvg = nan(nBins,1);
valMid = nan(nBins,1);
for ii = 1 : nBins
    valMid(ii) = mean(bins(ii:ii+1));
end


%Write range as variable output argument
if nargout > 2
   varargout{1} = rng; 
   
   if nargout > 3
       varargout{1} = valAvg; 
    end
end
    
if iscell(data)
   error('ePdf:cell','This function cannot handle cell input.'); 
end


%Sort data:
data(isnan(data)) = [];
data = sort(data(:));

if numel(data) == 0
   return
end
if all(isnan(data))
    error('ePdf:dataNan','All input data are nan. This is not allowed.');
end

%(Loop over bins)
for ii = 1 : nBins
    
    if ii < nBins
        indCurr = (find(data >= bins(ii), 1, 'first') : find(data < bins(ii+1), 1, 'last'));
    else
        indCurr = (find(data >= bins(ii), 1, 'first') : find(data <= bins(ii+1), 1, 'last'));
    end
    
    if ~isempty(indCurr)
        prob(ii) = numel(indCurr);
        valAvg(ii) = nanmean(data(indCurr));
    else
        if regexpbl(noVal, 'zero')
            prob(ii) = 0;
        elseif regexpbl(noVal, 'nan')
            prob(ii) = nan;
        end

        valAvg(ii) = nan;
    end
end


if regexpbl(typ, 'frac')
    prob = prob / nPtUse;
elseif regexpbl(typ, 'percent')
    prob = 100 * prob / nPtUse;
elseif ~regexpbl(typ, 'count')
    error('ePdf:typUnknown',['The empirical PDF type ' typ ' has not been programmed for.']);
end

% %Calculate values corresponding to PDF:
% cdfIn = linspace(0,1,nData);
% cdfTemp = linspace(0,1,nData*nBins);
% dataTemp = interp1(cdfIn, data, cdfTemp, 'pchip');
% val = nan(bins,1);
% cntr = 0;
% for ii = 1 : bins
%     val(ii) = nanmean(dataTemp(cntr+1:cntr+nData));
%     cntr = cntr + nPtsIn;
% end
% clear cntr

% This check causes error in cases when the range is a subset of the input data 
%Check that PDF sums to 1
% if round2(sum(pdf), 4) ~= 1
%    error('ePdf:missingData',['The PDF sums to ' num2str(sum(pdf)) ',  but it should equal 1.0.']); 
% end



