function date = nc_days_2_date(sData, timeField)


if isfield(sData, 'atttime')
    attTime = sData.atttime;
elseif isfield(sData, 'attTime')
    attTime = sData.attTime;
else
    error('nd_days2date:attTime','No field was found for time attributes');
end

[gcmRef, gcmUnits] = NC_time_units(attTime);
cal =  NC_cal(attTime);

if ~regexpbl(gcmUnits, 'day')
    error('nd_days2date:timeUnits',['The time units are  ' gcmUnits ', which has not been prorammed for']);
end

if ~isfield(sData, timeField)
    error('nd_days2date:timeField',['The time field requested is ' ...
        timeField ', which is not present in the data structure.']);
end

date = days_2_date_v2(sData.(timeField), gcmRef, cal);