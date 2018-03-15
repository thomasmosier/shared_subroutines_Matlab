function [dataOutCurr, wgtAvgIn, wgtAvgOut] = norm_grid(dataInCurr, areaIn, dataOutCurr, areaOut, strNorm)

warning('normGrid:obsolete','This function is obsolete. Use norm_grid_v2');

% dataInCurr = squeeze(sEraDy.(varDs)(indDyCurr,:,:));
% dataOutCurr = squeeze(dataOutCurr);
ptsNormIn = find(~isnan(dataInCurr));
ptsNormOut = find(~isnan(dataOutCurr));
wgtAvgIn  = (sum2d(dataInCurr(ptsNormIn).*areaIn(ptsNormIn))/sum2d(areaIn(ptsNormIn)));
wgtAvgOut = (sum2d(dataOutCurr(ptsNormOut).*areaOut(ptsNormOut))/sum2d(areaOut(ptsNormOut)));
if regexpbl(strNorm, 'mult')
    %Normalize output by input (enforce conservation of mass):
    dataOutCurr = squeeze(dataOutCurr) / (wgtAvgOut/wgtAvgIn);
elseif regexpbl(strNorm, 'add')
    %Normalize output by input (enforce conservation of mass):
    dataOutCurr = squeeze(dataOutCurr) - (wgtAvgOut - wgtAvgIn);
elseif regexpbl(strNorm, 'none')
%
else
    error('ERAInterp:unknownNorm',['The normalization method ' strNorm ' has not been prgorammed for.'])
end