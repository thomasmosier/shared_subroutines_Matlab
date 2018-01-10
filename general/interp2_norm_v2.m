function [dataOut, varargout] = interp2_norm_v2(lonIn, latIn, areaIn, dataIn, lonOut, latOut, areaOut, strIntrp, strNorm)

%Do Interpolation
dataIn = squeeze(dataIn);

if regexpbl(strIntrp,{'area','conserve','remap'},'and') || regexpbl(strNorm,{'area','conserve','remap'},'and')
    dataOut = area_conserve_remap(lonIn, latIn, dataIn, lonOut, latOut);
else
    if regexpbl(strIntrp,'pchip')
            warning('off', 'PCHIP_2D:NaN');
        dataOut = PCHIP_2D(lonIn, latIn, dataIn, lonOut, latOut);
            warning('on', 'PCHIP_2D:NaN');
    else
        dataOut = interp2(lonIn, latIn, dataIn, lonOut, latOut, strIntrp);
    end


    %%Normalize output:
    [dataOut, wgtIn, wgtOut] = norm_grid_v2(dataIn, areaIn, latIn, lonIn, dataOut, areaOut, latOut, lonOut, nan(1,2), nan(1,2), strNorm);
end

if nargout == 3
    varargout{1} = wgtIn;
    varargout{2} = wgtOut;
end


% if regexpbl(strNorm, {'mean', 'avg'})
% %     norm = (mean2d(dataOut)) / (mean2d(dataIn));
%     norm = (sum2d(areaOut.*dataOut)/sum2d(areaOut)) / (sum2d(areaIn.*dataIn)/sum2d(areaIn));
% elseif regexpbl(strNorm, {'sum', 'vol'})
%     norm = (sum2d(areaOut.*dataOut)/sum2d(areaOut)) / (sum2d(areaIn.*dataIn)/sum2d(areaIn));
% elseif regexpbl(strNorm, {'none'})
%     norm = 1;
% else
%     error('interpNorm:unknownType', ['The nromalization method ' ...
%         strNorm ' has not been programmed for.']);
% end

