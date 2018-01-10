function [data, lat, lon] = rm_nan_border(data, lat, lon)
%Function to remove data elements that are all nan

%Remove nan's from left column
while all(isnan(data(:,1)))
    data = data(:,2:end);
    lon = lon(:,2:end);
end

%Remove nan's from right column
while all(isnan(data(:,end)))
    data = data(:,1:end-1);
    lon = lon(:,1:end-1);
end

%Remove nan's from top row
while all(isnan(data(1,:)))
    data = data(2:end,:);
    lat = lat(2:end,:);
end

%Remove nan's from bottom row
while all(isnan(data(end,:)))
    data = data(1:end-1,:);
    lat = lat(1:end-1,:);
end