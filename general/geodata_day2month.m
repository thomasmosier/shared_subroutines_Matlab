function sData = geodata_day2month(sData, var, type)

varDate = 'date';

szIn = size(sData.(var));
newDates = unique(sData.(varDate)(:,1:2), 'rows');

temp = nan([numel(newDates(:,1)), szIn(2:3)], 'single');
for ii = 1 : numel(newDates(:,1))
    indAgg = find(ismember(sData.(varDate)(:,1:2), newDates(ii,:), 'rows') == 1);
    if strcmpi(type, 'sum')
        temp(ii,:,:) = squeeze(sum(sData.(var)(indAgg,:,:), 1));
    elseif strcmpi(type, 'avg')
        temp(ii,:,:) = squeeze(mean(sData.(var)(indAgg,:,:), 1));
    else
        error('geodataday2month:unknownType', ['The aggregation type ' type ' has not been pgorammed for.']);
    end
end
sData.(var) = temp;
sData.(varDate) = newDates;
sData = rmfield(sData, 'time');