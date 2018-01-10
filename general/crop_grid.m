function [data, lon, lat] = crop_grid(data, lon, lat, lonCrp, latCrp, fr, type)

%Find indices to crop to:
[indLonGp, indLatGp] = find_crop_ind(lon, lat, lonCrp, latCrp, fr, type);

%Crop
lon = lon(indLonGp(1) : indLonGp(end));
lat  = lat( indLatGp(1) : indLatGp(end));
if ndims(data) == 3
    data = data(:, indLatGp(1) : indLatGp(end), indLonGp(1): indLonGp(end));
elseif ismatrix(data)
    data = data(   indLatGp(1) : indLatGp(end), indLonGp(1): indLonGp(end));
end