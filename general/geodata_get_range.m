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

mn = double(min(dataTmp));
mx = double(max(dataTmp));
avg = double(mean(dataTmp));
sd = double(std(dataTmp));