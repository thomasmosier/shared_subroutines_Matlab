function sDataClim = geodata_clim(sData, varData, yrsClim)


%Find time-series elements to compare during bias correction:
%For historic model:
if isfield(sData, 'atttime')
    timeAtt = 'atttime';
elseif isfield(sData, 'attTime')
    timeAtt = 'attTime';
else
    error('eQMGeodata:noTimeAttMod','No time attribute field found in the historic modelled data.')
end
[refData, timeUnits] = NC_time_units(sData.(timeAtt));
if regexpbl(timeUnits, 'days')
    dateIn = days_2_date(sData.time, refData, 'gregorian');
else
    error('eQMGeodata:unknownTimeUnit','The time unit for the historic model is not known.')
end

indUse = find(ismember(dateIn(:,1), (yrsClim(1):yrsClim(end))','rows') > 0);

%Calculate output
sDataClim = sData;
sDataClim.data = squeeze(mean(sData.(varData)(indUse,:,:),1));
    sDataClim.time = NaN;
