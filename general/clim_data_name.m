function [nm, nmDisp] = clim_data_name(input)

% nm = 'unknown';
% nmDisp = 'unknown';

if regexpbl(input,{'world','clim'},'and')
    nm = 'WorldClim';
    nmDisp = 'WorldClim';
elseif regexpbl(input,'PRISM')
    nm = 'PRISM';
    nmDisp = 'PRISM';
elseif regexpbl(input,{'will','matsuura','udel','W&M'})
    nm = 'WM';
    nmDisp = 'Willmott and Matsuura';
elseif regexpbl(input,{'CRU'}) || regexpbl(input,{'climate','research','unit'},'and')
    nm = 'CRU';
    nmDisp = 'Climate Research Unit';
elseif regexpbl(input,'GPCC') || regexpbl(input,{'global','precipitation','climatology'},'and')
    nm = 'GPCC';
    nmDisp = 'Global Precipitation Climatology Centre';
elseif regexpbl(input,'era')
    nm = 'ERA';
    nmDisp = 'ERA';
elseif regexpbl(input,'srtm')
    nm = 'srtm';
    nmDisp = 'Shuttle Radar Topography Mission';
else
    [~, nm, ~] = fileparts(input);
    indHG = regexpi(nm, '_');
    if isempty(indHG)
        indHG = regexpi(nm, '\.');
    end

    if ~isempty(indHG)
        if regexpbl(nm,{'Amon','historic'}) && numel(indHG) > 1
            nm = nm(indHG(2)+1:indHG(end)-1);
            nmDisp = ['CMIP5 GCM ' nm];
        else
            nmDisp = nm;
            warning('main:input_TS_type', ['The program was not able to '...
                'detect the name of the current input: ' input]);
        end
    else
        nm = 'unknown';
        nmDisp = 'unknown';
    end
end