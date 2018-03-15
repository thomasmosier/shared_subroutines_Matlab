function sData = geodata_time_crop(sData, varOut, yrs, mnths)

varDate = 'date';
varTime = 'time';

if iscell(sData)
    for kk = 1 : numel(sData(:))
        indTimeUse = intersect( ...
            find(ismember(sData{kk}.(varDate)(:,2), mnths) == 1), ...
            find(sData{kk}.(varDate)(:,1) >= min(yrs) & sData{kk}.(varDate)(:,1) <= max(yrs)));

        sData{kk}.(varDate) = sData{kk}.(varDate)(indTimeUse,:);
        sData{kk}.(varOut) 	= sData{kk}.(varOut )(indTimeUse,:,:);
        if isfield(sData{kk}, varTime)
            sData{kk}.(varTime) = sData{kk}.(varTime)(indTimeUse);
        end
    end
    clear kk
elseif isstruct(sData)
    indTimeUse = intersect( ...
        find(ismember(sData.(varDate)(:,2), mnths) == 1), ...
        find(sData.(varDate)(:,1) >= min(yrs) & sData.(varDate)(:,1) <= max(yrs)));

    sData.(varDate) = sData.(varDate)(indTimeUse,:);
    sData.(varOut) 	= sData.(varOut )(indTimeUse,:,:);
    if isfield(sData, varTime)
        sData.(varTime) = sData.(varTime)(indTimeUse);
    end
else
    error('geodataTime:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.'])
end