function sData = geodata_spatial_crop(sData, varOut, latCrp, lonCrp, fr, type)

varLon = 'longitude';
varLat = 'latitude';

if isempty(latCrp) || all(isnan(latCrp))
    latCrp = [-1, 1]*10^6;
elseif isnan(latCrp(1))
    latCrp(1) = -10^6;
elseif isnan(latCrp(2))
   latCrp(2) = 10^6;
end

if isempty(lonCrp) || all(isnan(lonCrp))
    lonCrp = [-1, 1]*10^6;
elseif isnan(lonCrp(1))
    lonCrp(1) = -10^6;
elseif isnan(lonCrp(2))
   lonCrp(2) = 10^6;
end



if iscell(sData)
    for kk = 1 : numel(sData(:))
        [sData{kk}.(varOut), sData{kk}.(varLon), sData{kk}.(varLat)] ...
            = crop_grid(sData{kk}.(varOut), sData{kk}.(varLon), sData{kk}.(varLat), lonCrp, latCrp, fr, type);
    end
    clear kk
elseif isstruct(sData)
[sData.(varOut), sData.(varLon), sData.(varLat)] ...
            = crop_grid(sData.(varOut), sData.(varLon), sData.(varLat), lonCrp, latCrp, fr, type);
else
    error('geodataTime:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.'])
end