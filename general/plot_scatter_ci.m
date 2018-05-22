function plot_scatter_ci(sPlot, clLgd, varargin)

%Time-series are input as set of 3 arrays:
%(1) y-data (n by 1)
%(2) confidence intervals (n by 2)
%(3) x-data (n by 1)

if ischar(clLgd)
    clLgd = {clLgd};
end

nTs = numel(clLgd);

mrkrClr = [];
path = [];
yLab = [];
xLab = [];
ciType = 'err';
refLn = {};
mrkrSpec = [];
mrkrSize = [];
xLim = nan(2,1);
yLim = nan(2,1);
if numel(varargin(:)) > nTs*2
    for ii = nTs*2 + 1 : numel(varargin(:))
        if strcmpi(varargin{ii}, 'refline')
            refLn = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'color')
            mrkrClr = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'path')
           path = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'ylab')
            yLab = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'xlab')
            xLab = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'citype')
            ciType = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'mrkrspec')
            mrkrSpec = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'size')
            mrkrSize = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'xlim')
            xLim = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'ylim')
            yLim = varargin{ii+1};
        end
    end
end

if isempty(mrkrClr)
    mrkrClr = distinguishable_colors(nTs);
else
    if numel(mrkrClr(:,1)) ~= nTs 
        error('plotTsCi:clrNumber','The number of colors provided does not match the number of time-series detected.')
    end
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

%Plot error bars first:
hCI = nan(nTs, 1);
for ii = 1 : nTs
%     if ~isempty(varargin{2+3*(ii-1)})
%         if strcmpi(ciType, 'err')
%             ciCurr = varargin{2+3*(ii-1)} + repmat(varargin{1+3*(ii-1)}, [1,2]);
%         elseif strcmpi(ciType, 'abs')
%             ciCurr = varargin{2+3*(ii-1)};
%         else
%             error('plotTsCi:unknownCiType',['The confidence interval type ' ciType ' is unknwown.']);
%         end
%         
%         
%         hCI(ii) = ciplot(ciCurr(:,1), ciCurr(:,2), varargin{3+3*(ii-1)}, tsClr(ii,:));
%         
%         alpha(hCI(ii), 0.2);
%         
%         yMn = min(yMn, min(ciCurr(:,1)));
%         yMx = max(yMx, max(ciCurr(:,2)));
%     end
end

hCI(isnan(hCI)) = [];

if ~isempty(mrkrSpec)
    if ~iscell(mrkrSpec)
        mrkrSpec = {mrkrSpec};
    end
    
    if numel(mrkrSpec) == 1
        mrkrSpec = repmat(mrkrSpec, nTs, 1);
    elseif ~numel(mrkrSpec) == nTs
        error('plotTsCi:nLineSpec', ['The line specification has ' num2str(numel(mrkrSpec)) ' entries, but there are ' num2str(nTs) ' lines to plot.'])
    end
end

%Plot series:
hScat = nan(nTs,1);
for ii = 1 : nTs
    xCurr = varargin{1+2*(ii-1)};
    yCurr = varargin{2+2*(ii-1)};

    if ~isempty(mrkrSpec)
        hScat(ii) = scatter(xCurr, yCurr, mrkrSize, mrkrClr(ii,:), mrkrSpec{ii}, 'filled');
    else
        hScat(ii) = scatter(xCurr, yCurr, mrkrSize, mrkrClr(ii,:), 'filled');
    end
    
    yMn = min(yMn, min(yCurr));
    yMx = max(yMx, max(yCurr));
    
    xMn = min(xMn, min(xCurr));
    xMx = max(xMx, max(xCurr));
end


if ~isempty(refLn(:))
    hLine = nan(numel(refLn));

    for ii = 1 : numel(refLn)
        hLine(ii) = line(refLn{ii}(1:2), refLn{ii}(3:4));
        refGray = [0.5,0.5,0.5];

        set(hLine(ii),'LineWidth', sPlot.lnwd, ...
            'color', refGray, 'LineStyle','--');
    end
end

hold off

%Set x and y limits (if input)
if all(~isnan(xLim))
    xlim(xLim);
else
    xRng = xMx - xMn;
    if xRng  < 10
        xScl = 1.07;
    elseif xRng < 20
        xScl = 1.07;
    elseif xRng  < 30
        xScl = 1.09;
    elseif xRng  < 40
        xScl = 1.09;
    else
        xScl = 1.1;
    end
    xlim([0.97*xMn, xScl*xMx]);
end
if all(~isnan(yLim))
    ylim(yLim);
else
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
    ylim([0.97*yMn, yScl*yMx]);
end

%Check to see if there are legend entries for all inputs:
hLgd = [];
if numel(clLgd(:)) == numel(hScat(:))
    hTsLgd = hScat;
    for ii = numel(hScat(:)) : -1 : 1
        if isempty(clLgd{ii})
            clLgd(ii) = [];
            hTsLgd(ii) = [];
        end
    end
    
    if ~isempty(hTsLgd)
        hLgd = legend(hTsLgd, clLgd, 'Location','northwest');
    end
else
    warning('plotTsCi:diffNumTsLgd',['There are ' num2str(numel(hScat(:))) ...
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
set(hScat, ...
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
