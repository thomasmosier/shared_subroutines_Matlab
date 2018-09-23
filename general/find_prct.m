function [indOut, varargout] = find_prct(probIn, valIn, prcThresh, pdfTyp)


probNNan = flipud(probIn(~isnan(probIn)));
pdfCum = cumsum(probNNan);

if numel(pdfTyp, 'count')
    [indOut] = find(pdfCum < prcThresh*max(pdfCum)/100, 1, 'last');
elseif numel(pdfTyp, 'percent')
    [indOut] = find(pdfCum < 100 - prcThresh, 1, 'last');
elseif numel(pdfTyp, 'frac')
    [indOut] = find(pdfCum < 1 - prcThresh/100, 1, 'last');
end

indOut = (find(probNNan(indOut) == probIn) : numel(probIn));

if nargout > 1
    probOut = nan(size(probIn));
    probOut(indOut) = probIn(indOut);
    
    varargout{1} = probOut;
    if nargout > 2
        valOut = nan(size(valIn));
        valOut(indOut) = valIn(indOut);
    
        varargout{2} = valOut;
    end
end