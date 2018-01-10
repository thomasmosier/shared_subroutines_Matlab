function sSrtm = read_SRTM_tif(pathSrtm, varElev, lonDs, latDs, fr, cropType)


sSrtm = struct;

%METHOD USING 'geoimread' (From File Exchange)
[sSrtm.(varElev),sSrtm.longitude,sSrtm.latitude,~] = geoimread(pathSrtm);
    sSrtm.latitude = sSrtm.latitude(:);
    
% %METHOD USING MATLAB's BUILt-IN 'geotiffread' (contained in the mapping
% %toolbox)
% [srtmIn, srtmR] = geotiffread(pathSrtm);
% %Assign to structure array:
% 
% sSrtm.longitude = (srtmR.LongitudeLimits(1) + 0.5*srtmR.CellExtentInLongitude ...
%     : srtmR.CellExtentInLongitude ...
%     : srtmR.LongitudeLimits(2) - 0.5*srtmR.CellExtentInLongitude); 
% sSrtm.latitude = (srtmR.LatitudeLimits(2) - 0.5*srtmR.CellExtentInLatitude ...
%     : -srtmR.CellExtentInLatitude ...
%     : srtmR.LatitudeLimits(1) + 0.5*srtmR.CellExtentInLatitude)';
% sSrtm.(varElev) = srtmIn;

% %METHOD USING 'geotiff_read' (From File Exchange)
%DOES NOT WORK PROPERLY!!!
% I = geotiff_read_v2(pathSrtm);

%Crop:
sSrtm = crop_geostruct(sSrtm, varElev, lonDs, latDs, fr, cropType);

%No data values: 32767
valNoData = 32767;
sSrtm.(varElev)(sSrtm.(varElev) == valNoData) = nan;