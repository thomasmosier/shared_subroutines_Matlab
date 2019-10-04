function sData = read_geotif(pathData, varLd, lonLd, latLd, fr, cropType)


sData = struct;

%METHOD USING 'geoimread' (From File Exchange)
[sData.(varLd),sData.longitude,sData.latitude,I] = geoimread(pathData);
    sData.latitude = sData.latitude(:);
    
%If geotiff is in UTM, get UTM zone
if regexpbl(I.Projection, 'UTM')
    sData.utmZone = sprintf('%02d %s', I.Zone, I.Projection(end:end));
end

%Crop:
if all(~isnan(lonLd)) && all(~isnan(latLd))
    sData = crop_geostruct(sData, varLd, lonLd, latLd, fr, cropType);
end

%No data values: 32767
nanValTest = [min2d(sData.(varLd)), max2d(sData.(varLd))];
[~, indMx] = max(abs(nanValTest));
valNoData = nanValTest(indMx);

% disp(['The read_geotif function is setting the no-data value to ' num2str(valNoData) '.']);

sData.(varLd)(sData.(varLd) == valNoData) = nan;