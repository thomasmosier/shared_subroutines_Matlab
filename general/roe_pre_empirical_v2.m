function varargout = roe_pre_empirical_v2(vaporSat, demGradWind, W, mode, paramIn, varargin)

% vaporSat = squeeze(nanmean(vaporSat, 1));
% demGrad = squeeze(nanmean(demGrad, 1));
% dem = squeeze(nanmean(dem, 1));
%     preRef = squeeze(nanmean(preRef, 1));
    
szIn = size(vaporSat);

% Hm = 3; %Units = km
% W = neighbor_wgt_matrix(size(vaporSat), 'perp');
T = double([vaporSat(:), demGradWind(:).*vaporSat(:)]);
% T = [vaporSat(:), demGradWind(:).*vaporSat(:).*exp(-dem(:)/Hm)];
        
if regexpbl(mode, 'cal')
    preRef = varargin{1};

    parameters = combvec(paramIn{1}(:)', paramIn{2}(:)', paramIn{3}(:)', paramIn{4}(:)');
    
    nComb = numel(parameters(1,:));
    maeResidual = nan(nComb, 1);
    biasResidual = maeResidual;
    for xx = 1 : nComb
        preMod = reshape(T*parameters(1:2,xx) + W*T*parameters(3:4,xx), szIn);
        preMod(preMod < 0) = 0; 
        
        maeResidual(xx) = nanmean(abs(preMod(:) - preRef(:)));
        biasResidual(xx) = nanmean(preMod(:)) - nanmean(preRef(:));
    end
    clear xx

    varargout{1} = parameters';
    varargout{2} = maeResidual;
    varargout{3} = biasResidual;
elseif regexpbl(mode, 'val')
    if ~iscell(paramIn) || numel(paramIn{1}) ~= 1 
        error('roePreEmpirical:incorrectInput','The input variables do not have the correct length');
    end

    parameters = [paramIn{1}, paramIn{2}, paramIn{3}, paramIn{4}];
    
    preCov = reshape(T(:,1)*parameters(1) + W*T(:,1)*parameters(3), szIn);
    preOro = reshape(T(:,2)*parameters(2) + W*T(:,2)*parameters(4), szIn);

%     preMod = reshape(T*parameters(1:2) + W*T*parameters(3:4), szIn);
    
    varargout{1} = preCov;
    varargout{2} = preOro;
end

% %First round of optimize orographic precipitation relationship:
% alpha = linspace( -2,  1, nPts);
% beta = linspace( -2,  1, nPts);
% Hm  = linspace(0.1, 10, nPts);
% maeResidual = nan(size([numel(alpha), numel(beta), numel(Hm)]));
% biasResidual = maeResidual;
% for xx = 1 : numel(alpha)
%     for yy = 1 : numel(beta)
%         for zz = 1 : numel(Hm)
%             alphaCurr = alpha(xx);
%             betaCurr = beta(yy); %units = kg/m^3
%             HmCurr = Hm(zz); %units = km
%             totPrecip = 10^(alphaCurr)*sEraDy.(varVapSat) + 10^(betaCurr)*sEraDy.(varHussSat)...
%                 .*sEraDy.(varWindDot).*exp(eraDemTemp/HmCurr);
%             totPrecip(totPrecip < 0) = 0; 
%             maeResidual(xx,yy,zz) = nanmean(abs(totPrecip(:) - sEraDy.(varDs)(:)));
%             biasResidual(xx,yy,zz) = nanmean(totPrecip(:)) - nanmean(sEraDy.(varDs)(:));
%         end
%     end
% end
% clear xx yy zz