function plot_ts(pathOut, dataIn, datesIn, namesIn, yrsUse, tstep, type, var, units)

strVis = 'off';


%%Plot options:
%Set font size:
ftSz = 11;
%     ftSzT = ftSz + 1;
%     ftSzAx = ftSz - 1;
%Set line width:
lnWd = 2;
    lnWdA = lnWd - 0.5;
%     mrkrWd = lnWd;
strFont = 'Arial';
% %Use these colors instead?
% %From color Brewer):
% clrBrewer = [228,26,28; ...
%     55,126,184; ...
%     77,175,74; ...
%     152,78,163]/255;
%Dimensions of plot frame:
szFrameLg = [8.5,6]; %width, height
% szFrameWd = [18, 6]; %width, height
% szFrameSm = [5,3];
res = 300;

%Turn off image export warning:
warning('off', 'MATLAB:LargeImage');


%If data not cell, put in cell:
if isnumeric(dataIn) && isnumeric(datesIn)
    dataTemp = dataIn;
    datesTemp = datesIn;
    
    dataIn = cell(numel(dataIn(1,:)), 1);
    datesIn = dataIn;
    for ii = 1 : numel(datesIn(:))
        dataIn{ii} = dataTemp(:,ii);
        datesIn{ii} = datesTemp;
    end
end



nSeries = numel(dataIn(:));

%Make figures of annual, monthly, or daily values
lsYrs = (min(yrsUse):1:max(yrsUse));
nYrs = numel(lsYrs);
if regexpbl(tstep, {'ann','year'})
    indDateChck = 1;
    strXLab = 'Year';
    strYLab = 'Annual';
elseif regexpbl(tstep, {'month','mnth'})
    indDateChck = 2;
    strXLab = 'Month/Year';
    strYLab = 'Monthly';
elseif regexpbl(tstep, {'day','daily'})
    indDateChck = 3;
    strXLab = 'Day/Month/Year';
    strYLab = 'Daily';
else
    error('plotTs:unknownTimeStep',['The time step ' tStep ' has not been programmed for.']);
end


%Find all dates present:
datesTemp = nan(0, indDateChck);
for jj = 1 : nSeries
    datesTemp = unique(union(datesTemp, datesIn{jj}(:,1:indDateChck), 'rows'), 'rows');
end
datesTemp = datesTemp( datesTemp(:,1) >= min(yrsUse) & datesTemp(:,1) <= max(yrsUse), 1:indDateChck);

datesTemp = sortrows(datesTemp);

%Fill in dates to ensure continuous record:
yrsIn = (min(datesTemp(:,1)):max(datesTemp(:,1)))';
if regexpbl(tstep, {'ann','year'})
    datesPlot = yrsIn;
elseif regexpbl(tstep, {'month','mnth'})
    mnthsIn = unique(datesTemp(:,2));
    
    nMnths = numel(mnthsIn);
    datesPlot = nan(0, indDateChck);
    for ii = 1 : numel(yrsIn)
        datesPlot(end+1:end+numel(mnthsIn), :) = [repmat(yrsIn(ii), [nMnths, 1]), mnthsIn(:)];
    end
elseif regexpbl(tstep, {'day','daily'})
    mnthsIn = unique(datesTemp(:,2));
    
    nMnths = numel(mnthsIn);
    datesPlot = nan(0, indDateChck);
    for ii = 1 : numel(yrsIn)
        for jj = 1 : nMnths
            nDys = eomday(yrsIn(ii), mnthsIn(jj));
            datesPlot(end+1:end+nDys, :) = [repmat(yrsIn(ii), [nDys, 1]), repmat(mnthsIn(jj), [nDys, 1]), (1:nDys)'];
        end
    end
else
    error('plotTs:unknownTimeStep',['The time step ' tStep ' has not been programmed for.']);
end

nTsPts = numel(datesPlot(:,1));

dataPlot = nan(nTsPts, nSeries);
for ii = 1 : nTsPts
    for jj = 1 : nSeries
        indCurr = find(ismember(datesIn{jj}(:,1:indDateChck),datesPlot(ii,1:indDateChck), 'rows') == 1);
        
        if ~isempty(indCurr)
            if regexpbl(type, {'total','sum'})
                dataPlot(ii,jj) = sum(dataIn{jj}(indCurr));
            elseif regexpbl(type, {'avg','average','mean'})
                dataPlot(ii,jj) = mean(dataIn{jj}(indCurr));
            else
               error('plotTsAnn:type',['Type ' type ' is not known.']) 
            end
        end
    end
    
end

%Remove series that are entirely nan
for jj = nSeries : -1 : 1
    if sum(~isnan(dataPlot(:,jj)))/numel(dataPlot(:,jj)) < 0.03 %all(isnan(dataPlot(:,jj)))
       dataPlot(:,jj) = [];
       namesIn(jj) = [];
       nSeries = nSeries - 1;
    end
end

indXPlot = (1:nTsPts);

strDatesPlot = cell(nTsPts,1);
for ii = 1 : nTsPts
    if regexpbl(tstep, {'ann','year'})
        strDatesPlot{ii} = num2str(datesPlot(ii,1));
    elseif regexpbl(tstep, {'month','mnth'})
        strDatesPlot{ii} = [num2str(datesPlot(ii,2)) '/' num2str(datesPlot(ii,1))];
    elseif regexpbl(tstep, {'day','daily'})
        strDatesPlot{ii} = [num2str(datesPlot(ii,3)), '/', num2str(datesPlot(ii,2)) '/' num2str(datesPlot(ii,1))];
    else
        error('plotTs:unknownTimeStep',['The time step ' tStep ' has not been programmed for.']);
    end
end

%%Annual time-series
hFigTs = figure('Units','in','Position',[2 2 szFrameLg], ...
    'paperunits','in','paperposition',[2 2 szFrameLg], 'color', 'white', 'visible', strVis);

hold on
hPltTs = nan(nSeries,1);
for jj = 1 : nSeries
    hPltTs(jj) = plot(indXPlot(:), dataPlot(:,jj), 'linewidth', lnWd);
end
hold off

%Change X Tick Labels:
ax = gca;
ax.XTick = round(linspace(1,nTsPts, 10));
ax.XTickLabel = strDatesPlot(ax.XTick);

xlim([min(indXPlot), max(indXPlot)]);


hTsXLab = xlabel(strXLab);
hTsYLab = ylabel([strYLab ' ' var ' (' units ')']);

for ii = 1 : numel(namesIn(:))
    namesIn{ii} = strrep(namesIn{ii}, '_', ' ');
end
   
hTsLgd = legend(hPltTs, namesIn);

set([hTsLgd; hTsXLab; hTsYLab], ...
    'FontName'   , strFont, ...
    'FontSize'   , ftSz);
set(gca, ...
    'FontName'    , strFont, ...
    'FontSize'    , ftSz, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.02 .02] , ...
    'LineWidth'   , lnWdA     , ...
    'FontSize'    , ftSz );
%Export:
[foldOut, fileOut, ~] = fileparts(pathOut);

if ~exist(foldOut,'dir')
    mkdir(foldOut);
end

savefig(hFigTs, fullfile(foldOut, [fileOut  '.fig']));
%         export_fig([pathPlot '.eps'],'-painters');
export_fig(fullfile(foldOut, [fileOut  '.png']),'-painters', ['-r' num2str(res)]);
