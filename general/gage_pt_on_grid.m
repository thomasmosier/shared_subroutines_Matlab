function pt = gage_pt_on_grid(lonGage, latGage, lon, lat)

[~, gageRow] = min(abs(latGage-lat));
[~, gageCol] = min(abs(lonGage-lon));
if latGage > max(lat) + 0.5*abs(diff(lat(1:2))) || latGage < min(lat) - 0.5*abs(diff(lat(end-1:end))) 
    warning('CCHF_backbone:gagePtRow','The latitutde of the gauge point may be outside the area being modeled.');
elseif  lonGage < min(lon) - 0.5*abs(diff(lon(1:2))) || lonGage > max(lon) + 0.5*abs(diff(lon(end-1:end))) 
    warning('CCHF_backbone:gagePtCol','The longitude of the gauge point may be outside the area being modeled.');
end

pt = round(sub2ind([numel(lat), numel(lon)], gageRow, gageCol));