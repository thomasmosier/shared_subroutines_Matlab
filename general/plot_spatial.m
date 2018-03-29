function plot_spatial(data, sPlot, pathFig, clrLabel, varargin)

ci = [];
lat = [];
lon = [];
if ~isempty(varargin(:))
   for ii = 1 : numel(varargin(:)) 
       if strcmpi(varargin{ii}, 'ci')
          ci = varargin{ii+1};
       elseif strcmpi(varargin{ii}, 'lat')
           lat = varargin{ii+1};
       elseif strcmpi(varargin{ii}, 'lon')
           lon = varargin{ii+1};    
       end
   end
   clear ii
end

  
%Create figure
hFig = figure('Units','in','Position',[2 2 sPlot.sz], 'paperunits','in','paperposition',[2 2 sPlot.sz], 'color', 'white', 'visible', sPlot.vis);
hold on

if ~isempty(lat) && ~isempty(lon)
    hPlot = imagesc(lon, lat, data);
else
    hPlot = imagesc(data);
end

set(gca,'YDir','normal')
set(hPlot,'AlphaData',~isnan(data));
whitebg([0.6,0.6,0.6]);
mxPerClrBar = max2d(abs(data));
valCBar = [-mxPerClrBar, mxPerClrBar];
hClr = colorbar; 
if isfield(sPlot, clrmap)
    colormap(sPlot.clrmap); 
end
caxis(valCBar);
set(hFig, 'color', 'white');
%set(hPDiff, 'EdgeColor', 'none');
%Add hatches:
if ~isempty(ci)
    if isempty(lat) || isempty(lon)
       error('plotSpatial:nocrd','Coordinates must be input when the optional confidence region is being used.'); 
    end
    
    ci(ci > 0) = 1;
    ci(ci < 0) = 0;
    % plot hatching region:
    indLwCf = find(ci == 0);
    [rLwCf, cLwCf] = ind2sub(size(ci), indLwCf);
    hLwCf = nan(numel(indLwCf), 1);
    edgLon = box_edg(lon);
    edgLat = box_edg(lat);
    for kk = 1 : numel(indLwCf)
        xCurr = [edgLon(cLwCf(kk)), edgLon(cLwCf(kk)+1), edgLon(cLwCf(kk)+1), edgLon(cLwCf(kk)  )];
        yCurr = [edgLat(rLwCf(kk)), edgLat(rLwCf(kk)  ), edgLat(rLwCf(kk)+1), edgLat(rLwCf(kk)+1)];
        hLwCf(kk) = patch('XData',xCurr, 'YData', yCurr);
        set(hLwCf(kk),'linestyle','none');
        hatchfill(hLwCf(kk));
    end
end
hold off; 
xlim([edgLon(indLonBnd(1)), edgLon(indLonBnd(end)+1)]);
ylim([edgLat(indLatBnd(end)+1), edgLat(indLatBnd(1))]);
hXLab = xlabel('Degrees East');
hYLab = ylabel('Degrees North');
hCLab = ylabel(hClr, clrLabel);

%             set(hTtl, ...
%                 'FontSize'   , ftSzT, ...
%                 'Color', 'black', ...
%                 'FontName'   , strFont);
set([hXLab, hYLab, hCLab], ...
    'FontSize'   , sPlot.ftsz, ...
    'Color', 'black', ...
    'FontName'   , sPlot.font);
set(gca, ...
    'Box'         , 'off', ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'YMinorTick'  , 'on'      , ...
    'FontName'   , sPlot.font, ...
    'fontSize'    , sPlot.axfntsz, ...
    'LineWidth'   , sPlot.lnwd);


if exist([pathFig, '.fig'], 'file')
    delete([pathFig, '.fig'])
end
if exist([pathFig, '.png'], 'file')
    delete([pathFig, '.png'])
end
if exist([pathFig, '.eps'], 'file')
    delete([pathFig, '.eps'])
end
if exist([pathFig, '.tiff'], 'file')
    delete([pathFig, '.tiff'])
end
savefig(hFig, [pathFig '.fig']);
export_fig([pathFig '.eps'],'-painters', strFigRes);
export_fig([pathFig '.png'],'-painters', strFigRes);
export_fig([pathFig '.tiff'],'-painters', strFigRes);