function sData = geodata_ocean2nan(sData, varOut, sElev, varElev)


if iscell(sData) && iscell(sElev)
    for kk = 1 : numel(sData(:))
        indOcean = find(sElev{kk}.(varElev) == 0);
        [rOcean, cOcean] = ind2sub(size(sElev{kk}.(varElev)), indOcean);
    
        for ll = 1 : numel(indOcean)
            sData{kk}.(varOut)(:,rOcean(ll),cOcean(ll)) = nan;
        end
        clear ll
    end
    clear kk
elseif iscell(sData) && isstruct(sElev)
    %Set all values where elevation = 0 to nan (ocean)
    indOcean = find(sElev.(varElev) == 0);
    [rOcean, cOcean] = ind2sub(size(sElev.(varElev)), indOcean);
    % indNocea = find(sElev.(varElev) ~= 0);
    % [rNocea, cNocea] = ind2sub(size(sElev.(varElev)), indNocea);

    for kk = 1 : numel(sData(:))
        for ll = 1 : numel(indOcean)
            sData{kk}.(varOut)(:,rOcean(ll),cOcean(ll)) = nan;
        end
        clear ll
    end
    clear kk
elseif isstruct(sData) && isstruct(sElev)
    %Set all values where elevation = 0 to nan (ocean)
    indOcean = find(sElev.(varElev) == 0);
    [rOcean, cOcean] = ind2sub(size(sElev.(varElev)), indOcean);
    % indNocea = find(sElev.(varElev) ~= 0);
    % [rNocea, cNocea] = ind2sub(size(sElev.(varElev)), indNocea);

    for ll = 1 : numel(indOcean)
        sData.(varOut)(:,rOcean(ll),cOcean(ll)) = nan;
    end
    clear ll
else
    error('geodataOceanNan:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.']);
end