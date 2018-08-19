function [mn, mx, avg, sd] = geodata_get_range(sData, var, typ)


if iscell(sData)
    nCell = numel(sData(:));
    
    dataTmp = [];
    for ii = 1 : nCell
        dataTmp = [dataTmp; sData{ii}.(var)(:)];
    end
else
    dataTmp = sData.(var)(:);
end

if regexpbl(typ, 'nan')
    dataTmp(isnan(dataTmp)) = [];
elseif ~regexpbl(typ, 'all')
    error('geodataGetRange:typeUnknown', ['The type ' typ ' is unknown. It should either be nan or all.']);
end

mn = min(dataTmp);
mx = max(dataTmp);
avg = mean(dataTmp);
sd = std(dataTmp);