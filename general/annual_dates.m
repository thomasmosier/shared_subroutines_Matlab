function [indStrt, indEnd] = annual_dates(dateVec, dateStrtIn)



if numel(dateStrtIn) == 3 
    dateStrt = dateStrtIn(2:3);
elseif numel(dateStrtIn) == 2
    dateStrt = dateStrtIn;
elseif numel(dateStrtIn) ~= 2
    error('annualDates:wrongsizeDateStart', ['The datestart vector has ' num2str(numel(dateStrtIn)) '. 2 or 3 is expected.'])
end

if numel(dateVec(1,:)) ~= 3
    error('annualDates:wrongsizedateVec', ['The datevec vector has ' num2str(numel(dateStrtIn)) '. 2 or 3 is expected.'])
end

%Find unique years:
indStrt = find(ismember(dateVec(:,2:3), dateStrt, 'rows'));
    indStrt = indStrt(:); %Ensure column vector
if ~isempty(indStrt)
    indEnd = [indStrt(2:end) - 1; numel(dateVec(:,1))];
else
    warning('hydroStats:waterYearNoWY', ['Hydrologic statistics for ' ...
    sMeta.region{kk} ' are being skipped because no water years have been found.']);
end

if numel(dateStrtIn) == 3 
    indStrt = indStrt(1);
    indEnd = indEnd(1);
end


if numel(indEnd) > 1
    yrTest = nanmean(indEnd(1:end-1) - indStrt(1:end-1));
    if indEnd(end) - indStrt(end) < yrTest - 2 || indEnd(end) - indStrt(end) > yrTest + 2
        indStrt = indStrt(1:end-1);
        indEnd = indEnd(1:end-1);
    end
end