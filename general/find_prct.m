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
    varargout{1} = probIn(indOut);
    if nargout > 2
        varargout{2} = valIn(indOut);
    end
end