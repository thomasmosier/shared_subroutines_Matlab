function [dateIn, dataIn] = stn_csv_find_read_crop(foldHist, var, stn, yrs, mnths)

pathHistTmx = fullfile(foldHist, char(extract_field(dir(fullfile(foldHist, [var '_*.csv'])), 'name')));
[dateIn, dataIn] = stn_csv_read(pathHistTmx, stn);

indUse = find(dateIn(:,1) <= max(yrs) & dateIn(:,1) >= min(yrs) & ismember(dateIn(:,2), mnths) == 1);
dateIn = dateIn(indUse,:);
dataIn = dataIn(indUse,:);