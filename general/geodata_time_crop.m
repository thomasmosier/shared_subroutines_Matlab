function [sData, varargout] = geodata_time_crop(sData, varOut, yrs, mnths)

varDate = 'date';
varTime = 'time';

if isempty(yrs) || all(isnan(yrs))
    yrs = [-1, 1]*10^6;
elseif isnan(yrs(1))
    yrs(1) = -10^6;
elseif isnan(yrs(2))
   yrs(2) = 10^6;
end

if isempty(mnths) || all(isnan(mnths))
    mnths = (1:12);
end



if iscell(sData)
    indTimeUse = cell(numel(numel(sData(:))), 1);
    for kk = 1 : numel(sData(:))
        nDim = ndims(sData{kk}.(varOut));
        if numel(sData{kk}.(varDate)(1,:)) > 1
            indTimeUse{kk} = intersect( ...
                find(ismember(sData{kk}.(varDate)(:,2), mnths) == 1), ...
                find(sData{kk}.(varDate)(:,1) >= min(yrs) & sData{kk}.(varDate)(:,1) <= max(yrs)));
        elseif numel(sData{kk}.(varDate)(1,:)) == 1
            indTimeUse{kk} = find(sData{kk}.(varDate)(:,1) >= min(yrs) & sData{kk}.(varDate)(:,1) <= max(yrs));
        else
            error('geodataTimeCrop:date', ['The date vector has precision ' num2str(numel(sData{kk}.(varDate)(1,:))) ', which is not expected.'])
        end

        sData{kk}.(varDate) = sData{kk}.(varDate)(indTimeUse{kk},:);
        if nDim == 3
            sData{kk}.(varOut) 	= sData{kk}.(varOut )(indTimeUse{kk},:,:);
        else
            sData{kk}.(varOut) 	= sData{kk}.(varOut )(indTimeUse{kk},:);
        end
        if isfield(sData{kk}, varTime)
            sData{kk}.(varTime) = sData{kk}.(varTime)(indTimeUse{kk});
        end
    end
    clear kk
elseif isstruct(sData)
    if numel(sData.(varDate)(1,:)) > 1
        indTimeUse = intersect( ...
            find(ismember(sData.(varDate)(:,2), mnths) == 1), ...
            find(sData.(varDate)(:,1) >= min(yrs) & sData.(varDate)(:,1) <= max(yrs)));
    elseif numel(sData.(varDate)(1,:)) == 1
        indTimeUse = find(sData.(varDate)(:,1) >= min(yrs) & sData.(varDate)(:,1) <= max(yrs));
    else
        error('geodataTimeCrop:date', ['The date vector has precision ' num2str(numel(sData.(varDate)(1,:))) ', which is not expected.'])
    end

    sData.(varDate) = sData.(varDate)(indTimeUse,:);
    if ndims(sData.(varOut)) == 3
        sData.(varOut) 	= sData.(varOut )(indTimeUse,:,:);
    else
        sData.(varOut) 	= sData.(varOut )(indTimeUse,:);
    end
    if isfield(sData, varTime)
        sData.(varTime) = sData.(varTime)(indTimeUse);
    end
else
    error('geodataTime:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.'])
end

if nargout > 1
    varargout{1} = indTimeUse;
end
