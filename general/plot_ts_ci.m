function plot_ts_ci(sPlot, clLgd, varargin)

%Time-series are input as set of 3 arrays:
%(1) y-data (n by 1)
%(2) confidence intervals (n by 2)
%(3) x-data (n by 1)

if ischar(clLgd)
    clLgd = {clLgd};
end

nTs = numel(clLgd);

refLn = [];
tsClr = [];
path = [];
yLab = [];
xLab = [];
ciType = 'err';
lineSpec = [];
MarkerFaceColor = [];
MarkerEdgeColor = [];
xLm = nan(1,2);
if numel(varargin(:)) > nTs*3
    for ii = nTs*3 + 1 : numel(varargin(:))
        if strcmpi(varargin{ii}, 'refline')
            refLn = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'color')
            tsClr = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'path')
           path = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'ylab')
            yLab = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'xlab')
            xLab = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'citype')
            ciType = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'linespec')
            lineSpec = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'MarkerFaceColor')
            MarkerFaceColor = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'MarkerEdgeColor')
            MarkerEdgeColor = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'xlim')
            xLm = varargin{ii+1};
        end
    end
end

if isempty(tsClr)
    tsClr = distinguishable_colors(nTs);
else
    if numel(tsClr(:,1)) ~= nTs 
        error('plotTsCi:clrNumber','The number of colors provided does not match the number of time-series detected.')
    end
end

xMn = 10^6;
xMx = -10^6;
yMn = 10^6;
yMx = -10^6;

%Create figure
hFig = figure('Units','in','Position',[2 2 sPlot.sz], 'paperunits','in','paperposition',[2 2 sPlot.sz], 'color', 'white', 'visible', sPlot.vis, 'AlphaMap', 1);
% set(0,'defaultfigurecolor',[1 1 1])
% whitebg([1,1,1]);
set(gcf,'color','w');
hold on

%Plot confidence intervals first:
hCI = nan(nTs, 1);
for ii = 1 : nTs
    if ~isempty(varargin{2+3*(ii-1)})
        if strcmpi(ciType, 'err')
            ciCurr = varargin{2+3*(ii-1)} + repmat(varargin{1+3*(ii-1)}, [1,2]);
        elseif strcmpi(ciType, 'abs')
            ciCurr = varargin{2+3*(ii-1)};
        else
            error('plotTsCi:unknownCiType',['The confidence interval type ' ciType ' is unknwown.']);
        end
        
        
        hCI(ii) = ciplot(ciCurr(:,1), ciCurr(:,2), varargin{3+3*(ii-1)}, tsClr(ii,:));
        
        alpha(hCI(ii), 0.2);
        
        yMn = min(yMn, min(ciCurr(:,1)));
        yMx = max(yMx, max(ciCurr(:,2)));
    end
end

hCI(isnan(hCI)) = [];

if ~isempty(lineSpec)
    if ~iscell(lineSpec)
        lineSpec = {lineSpec};
    end
    
    if numel(lineSpec) == 1
        lineSpec = repmat(lineSpec, nTs, 1);
    elseif ~numel(lineSpec) == nTs
        error('plotTsCi:nLineSpec', ['The line specification has ' num2str(numel(lineSpec)) ' entries, but there are ' num2str(nTs) ' lines to plot.'])
    end
end

%Plot mean projection lines:
hTs = nan(nTs,1);
for ii = 1 : nTs
    xCurr = varargin{3+3*(ii-1)};
    tsCurr = varargin{1+3*(ii-1)};
    
    if all(~isnan(xLm))
        indUse = find(xCurr >= min(xLm) & xCurr <= max(xLm));
        
        xCurr = xCurr(indUse);
        tsCurr = tsCurr(indUse);
    end
    
    if ~isempty(lineSpec)
        hTs(ii) = plot(xCurr, tsCurr, lineSpec{ii}, 'color', tsClr(ii,:), 'LineWidth', sPlot.lnwd);
    else
        hTs(ii) = plot(xCurr, tsCurr,'color', tsClr(ii,:), 'LineWidth', sPlot.lnwd);
    end
    
    if ~isempty(MarkerFaceColor)
        if strcmpi(MarkerFaceColor, 'filled')
            set(hTs(ii), 'MarkerFaceColor', tsClr(ii,:));
        else
            set(hTs(ii), 'MarkerFaceColor', MarkerFaceColor);
        end
    end
    
    if ~isempty(MarkerEdgeColor)
        if isprop(hTs(ii), 'MarkerEdgeColor')
            if strcmpi(MarkerEdgeColor, 'same')
                set(hTs(ii), 'MarkerEdgeColor', tsClr(ii,:));
            else
                set(hTs(ii), 'MarkerEdgeColor', MarkerEdgeColor);
            end
        end
    end
    
    
    yMn = min(yMn, min(tsCurr));
    yMx = max(yMx, max(tsCurr));
    
    xMn = min(xMn, min(xCurr));
    xMx = max(xMx, max(xCurr));
end

if ~isempty(refLn)
    hLine = line(refLn(1:2), refLn(3:4));
    refGray = [0.5,0.5,0.5];
    
    set(hLine,'LineWidth', sPlot.lnwd, ...
        'color', refGray, 'LineStyle','--');
end

hold off

%Check to see if there are legend entries for all inputs:
hLgd = [];
if numel(clLgd(:)) == numel(hTs(:))
    hTsLgd = hTs;
    for ii = numel(hTs(:)) : -1 : 1
        if isempty(clLgd{ii})
            clLgd(ii) = [];
            hTsLgd(ii) = [];
        end
    end
    
    if ~isempty(hTsLgd)
        hLgd = legend(hTsLgd, clLgd, 'Location','northwest');
    end
else
    warning('plotTsCi:diffNumTsLgd',['There are ' num2str(numel(hTs(:))) ...
        ' time-series and ' num2str(numel(clLgd(:))) ' legend entries. '...
        'A legend is not being produced because the numbers do not match.']);
end

if isfield(sPlot, 'ylabel')
    hYLab = ylabel(sPlot.ylabel);

    set(hYLab, ...
    'FontName'   , sPlot.font, ...
    'FontSize'   , ftSzT, ...
    'Color', 'black');
end
if isfield(sPlot, 'xlabel')
    hXLab = xlabel(sPlot.xlabel);

    set(hXLab, ...
    'FontName'   , sPlot.font, ...
    'FontSize'   , ftSzT, ...
    'Color', 'black');
end

yRng = yMx - yMn;
if yRng  < 10
    yScl = 1.07;
elseif yRng < 20
    yScl = 1.07;
elseif yRng  < 30
    yScl = 1.09;
elseif yRng  < 40
    yScl = 1.09;
else
    yScl = 1.1;
end

if xMx > xMn
    xlim([xMn, xMx]);
elseif isnan(xMn) || isnan(xMx)
    return
end
if yMx > yMn
    ylim([0.97*yMn, yScl*yMx]);
elseif isnan(yMn) || isnan(yMx)
    return
end


if ~isempty(yLab)
    hYLab = ylabel(yLab);
else
    hYLab = [];
end
if ~isempty(xLab)
    hXLab = xlabel(xLab);
else
    hXLab = [];
end

% for kk = 1 : nSce
%     if kk == 1
%         clrLn = custGrn;
%     elseif kk == 2
%         clrLn = custPrp;
%     end
%     set(hCntryProj(kk),'LineWidth',lnWd, ...
%         'Color',clrLn, ...
%         'linestyle','--');
%     set(hCntryProj(kk),'LineWidth',lnWd, ...
%         'Color',clrLn, ...
%         'linestyle','--');
% end
%Set figure properties:
set(hTs, ...
    'LineWidth', sPlot.lnwd);
if ~isempty(hCI)
    set(hCI,...
        'linestyle', 'none');
end
if ~isempty(hLgd)
    set(hLgd, ...
        'FontSize'   , sPlot.ftsz, ...
        'LineWidth', sPlot.axlnwd);
end
if ~isempty(xLab)
    set(hXLab, ...
        'FontSize'   , sPlot.ftsz, ...
        'LineWidth', sPlot.axlnwd);
end
if ~isempty(yLab)
    set(hYLab, ...
        'FontSize'   , sPlot.ftsz, ...
        'LineWidth', sPlot.axlnwd);
end
set(gca, ...
    'Box'         , 'off', ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'YMinorTick'  , 'on'      , ...
    'fontSize'    , sPlot.axfntsz, ...
    'LineWidth'   , sPlot.axlnwd, ...
    'FontName'   , sPlot.font);

if ~isempty(path)
    if exist([path, '.fig'], 'file')
        delete([path, '.fig'])
    end
    if exist([path, '.png'], 'file')
        delete([path, '.png'])
    end
    if exist([path, '.eps'], 'file')
        delete([path, '.eps'])
    end
    if exist([path, '.tiff'], 'file')
        delete([path, '.tiff'])
    end
    savefig(hFig, [path '.fig']);
    export_fig([path '.eps'],'-painters', sPlot.res);
    export_fig([path '.png'],'-painters', sPlot.res);
    export_fig([path '.tiff'],'-painters', sPlot.res);
end
