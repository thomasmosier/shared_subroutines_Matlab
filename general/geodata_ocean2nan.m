function sData = geodata_ocean2nan(sData, varOut, sElev, varElev)


%Set all values where elevation = 0 to nan (ocean)
indOcean = find(sElev.(varElev) == 0);
[rOcean, cOcean] = ind2sub(size(sElev.(varElev)), indOcean);
% indNocea = find(sElev.(varElev) ~= 0);
% [rNocea, cNocea] = ind2sub(size(sElev.(varElev)), indNocea);

if iscell(sData)
    for kk = 1 : numel(sData(:))
        for ll = 1 : numel(indOcean)
            sData{kk}.(varOut)(:,rOcean(ll),cOcean(ll)) = nan;
        end
        clear ll
    end
    clear kk
elseif isstruct(sData)
    for ll = 1 : numel(indOcean)
        sData.(varOut)(:,rOcean(ll),cOcean(ll)) = nan;
    end
    clear ll
else
    error('geodataOceanNan:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.']);
end