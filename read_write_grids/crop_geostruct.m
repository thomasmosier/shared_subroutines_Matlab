function sData = crop_geostruct(sData, varData, lonCrp, latCrp, fr, type)

if isfield(sData, 'lon')
    varLon = 'lon';
elseif isfield(sData, 'longitude')
    varLon = 'longitude';
else
    error('cropgeostruct:lonNotFound','No longitude variable was found.');
end

if isfield(sData, 'lat')
    varLat = 'lat';
elseif isfield(sData, 'latitude')
    varLat = 'latitude';
else
    error('cropgeostruct:latNotFound','No longitude variable was found.');
end

%Crop data:
[sData.(varData), sData.(varLon), sData.(varLat)] ...
    = crop_grid(sData.(varData), sData.(varLon), sData.(varLat), ...
    lonCrp, latCrp, fr, type);
