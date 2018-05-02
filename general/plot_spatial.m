function plot_spatial(data, sPlot, clrLabel, varargin)

ci = [];
lat = [];
lon = [];
path = [];
box = [];
latArr    = cell(0,1); 
lonArr    = cell(0,1);
magLatArr = cell(0,1);
magLonArr = cell(0,1);
clrArr    = cell(0,1);
if ~isempty(varargin(:))
   for ii = 1 : numel(varargin(:)) 
       if strcmpi(varargin{ii}, 'ci')
          ci = varargin{ii+1};
       elseif strcmpi(varargin{ii}, 'lat')
           lat = varargin{ii+1};
       elseif strcmpi(varargin{ii}, 'lon')
           lon = varargin{ii+1}; 
       elseif strcmpi(varargin{ii}, 'path')
           path = varargin{ii+1}; 
       elseif strcmpi(varargin{ii}, 'box')
           box = varargin{ii+1}; %[lonW, lonE, latS, latN]
       elseif strcmpi(varargin{ii}, 'arrows')
           latArr{end+1,1} = varargin{ii+1}; 
           lonArr{end+1,1} = varargin{ii+2};
           magLatArr{end+1,1} = varargin{ii+3};
           magLonArr{end+1,1} = varargin{ii+4};
           clrArr{end+1} = varargin{ii+5};
       end
   end
   clear ii
end

  
%Create figure
hFig = figure('Units','in','Position',[2 2 sPlot.sz], 'paperunits','in','paperposition',[2 2 sPlot.sz], 'color', 'white', 'visible', sPlot.vis);
hold on

if ~isempty(lat) && ~isempty(lon)
    hPlot = imagesc(lon, lat, data);
    
    if isfield(sPlot, 'lonlimit')
        indLonBnd = [find(lon >= min(sPlot.lonlimit), 1, 'first'), find(lon <= max(sPlot.lonlimit), 1, 'last')];
    end
    if isfield(sPlot, 'latlimit')
        indLatBnd = [find(lat <= max(sPlot.latlimit), 1, 'first'), find(lat >= min(sPlot.latlimit), 1, 'last')];
    end
    
    edgLon = box_edg(lon);
    edgLat = box_edg(lat);
else
    hPlot = imagesc(data);
    indLonBnd = [];
    indLatBnd = [];
end

set(gca,'YDir','normal')
set(hPlot,'AlphaData',~isnan(data));
whitebg(hFig, [0.6,0.6,0.6]);
%Color bar properties:
hClr = colorbar; 
if isfield(sPlot, 'clrmap')
    colormap(sPlot.clrmap); 
end
if isfield(sPlot, 'cntrzero') && sPlot.cntrzero == 1
    if isfield(sPlot, 'clrrng') && ~isempty(sPlot.clrrng) && all(~isnan(sPlot.clrrng))
        mxPerClrBar = max(sPlot.clrrng);
    else
        mxPerClrBar = max2d(abs(data));
    end
    valCBar = [-mxPerClrBar, mxPerClrBar];
    caxis(valCBar);
elseif isfield(sPlot, 'clrrng') && ~isempty(sPlot.clrrng) && all(~isnan(sPlot.clrrng))
    caxis(sPlot.clrrng);
end

set(hFig, 'color', 'white');
%set(hPDiff, 'EdgeColor', 'none');
%Add hatches:
if ~isempty(ci)
    if isempty(lat) || isempty(lon)
       error('plotSpatial:nocrdci','Coordinates must be input when the optional confidence region is being used.'); 
    end
    
    ci(ci > 0) = 1;
    ci(ci < 0) = 0;
    % plot hatching region:
    indLwCf = find(ci == 0);
    [rLwCf, cLwCf] = ind2sub(size(ci), indLwCf);
    hLwCf = nan(numel(indLwCf), 1);
    for kk = 1 : numel(indLwCf)
        xCurr = [edgLon(cLwCf(kk)), edgLon(cLwCf(kk)+1), edgLon(cLwCf(kk)+1), edgLon(cLwCf(kk)  )];
        yCurr = [edgLat(rLwCf(kk)), edgLat(rLwCf(kk)  ), edgLat(rLwCf(kk)+1), edgLat(rLwCf(kk)+1)];
        hLwCf(kk) = patch('XData',xCurr, 'YData', yCurr);
        set(hLwCf(kk),'linestyle','none');
        hatchfill(hLwCf(kk));
    end
end

%Plot arrows (if requested)
if ~isempty(magLonArr) && ~isempty(magLatArr) && ~isempty(lonArr) && ~isempty(latArr)
    nQuiv = numel(magLonArr(:));
    for ii = 1 : nQuiv
        quiver(lonArr{ii}, latArr{ii}, squeeze(magLonArr{ii}), squeeze(magLatArr{ii}), 'color',  clrArr{ii});
    end
end


if ~isempty(box)
    if ~isempty(lat) || ~isempty(lon)
        boxDef = [box(1), box(3), box(2) - box(1), box(4) - box(3)]; %[xLeft, ySouth, width, height]
        hRect = rectangle('Position',boxDef);
        set(hRect, 'LineWidth', sPlot.lnwd);
    else
        error('plotSpatial:nocrdBox','Coordinates must be input when the optional box argument is being used.'); 
    end
end

if ~isempty(indLonBnd)
    xlim([edgLon(indLonBnd(1)), edgLon(indLonBnd(end)+1)]);
end
if ~isempty(indLatBnd)
    ylim([edgLat(indLatBnd(end)+1), edgLat(indLatBnd(1))]);
end
hold off;
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