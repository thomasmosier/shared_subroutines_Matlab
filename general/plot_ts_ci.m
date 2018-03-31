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
if numel(varargin(:)) > nTs*3
    for ii = nTs*2 + 1 : numel(varargin(:))
        if strcmpi(varargin{ii}, 'line')
            refLn = varargin{ii}+1;
        elseif strcmpi(varargin{ii}, 'color')
            tsClr = varargin{ii}+1;
       elseif strcmpi(varargin{ii}, 'path')
           path = varargin{ii+1};
        end
    end
end

if isempty(tsClr)
    tsClr = distinguishable_colors(nTs);
end

xMn = 10^6;
xMx = -10^6;
yMn = 10^6;
yMx = -10^6;

%Create figure
hFig = figure('Units','in','Position',[2 2 sPlot.sz], 'paperunits','in','paperposition',[2 2 sPlot.sz], 'color', 'white', 'visible', sPlot.vis, 'AlphaMap', 1);
set(0,'defaultfigurecolor',[1 1 1])
whitebg([1,1,1])
hold on

%Plot confidence intervals first:
hCI = nan(nTs, 1);
for ii = 1 : nTs
    ciCurr = varargin{2+3*(ii-1)} + repmat(varargin{1+3*(ii-1)}, [1,2]);
    hCI(ii) = ciplot(ciCurr(:,1), ciCurr(:,2), varargin{3+3*(ii-1)}, tsClr(ii,:));
    alpha(hCI(ii), 0.2);
    
    yMn = min(yMn, min(ciCurr(:,1)));
    yMx = max(yMx, max(ciCurr(:,2)));
    xMn = min(xMn, min(varargin{3+3*(ii-1)}(:)));
    xMx = max(xMx, max(varargin{3+3*(ii-1)}(:)));
end

%Plot mean projection lines:
hTs = nan(nTs,1);
for ii = 1 : nTs
    hTs(ii) = plot(varargin{3+3*(ii-1)}, varargin{1+3*(ii-1)},'color', tsClr(ii,:), 'LineWidth', sPlot.lnwd);
end

if ~isempty(refLn)
    hLine = line(refLn(1:2), refLn(2:4));
    refGray = [0.5,0.5,0.5];
    
    set(hLine,'LineWidth', sPlot.lnwd, ...
        'color', refGray, 'LineStyle','--');
end

hold off

hLgd = legend(hTs, clLgd, 'Location','northwest');

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
xlim([xMn, xMx]);
ylim([0.97*yMn, yScl*yMx]);


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
set(hCI,...
    'linestyle', 'none');
set(hLgd, ...
    'FontSize'   , sPlot.ftsz, ...
    'LineWidth', sPlot.axlnwd);
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
