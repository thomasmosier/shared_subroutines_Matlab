function sData = read_worldclim_tif(pathData, varLd, mnthsLd, lonDs, latDs, fr, cropType)

varLat = 'latitude';
varLon = 'longitude';

sData = struct;

%METHOD USING 'geoimread' (From File Exchange)
%Check if WordClim path includes file or is general directory
[root,file,extTemp] = fileparts(pathData);
%If includes file, remove so that it searches directory 
%(necessary because there are both hdr and data files):
if ~isempty(extTemp)
    pathWC = root;
else
    pathWC = fullfile(root,file);
end


if all(~isnan(mnthsLd))
    if regexpbl(varLd,{'tmp','tas','tmean'})   %Necessary because WorldClim uses 'tmean' instead of 'tmn'.
        wcVar = 'tavg';
    elseif regexpbl(varLd,{'tmn','tasmin','tmin'})
        wcVar = 'tmin';
    elseif regexpbl(varLd,{'tmx','tasmax','tmax'})
        wcVar = 'tmax';
    elseif regexpbl(varLd,'pre') || strcmpi(varLd,'pr')
        wcVar = 'prec';
    elseif regexpbl(varLd,{'alt','elev','DEM'})
        wcVar = 'alt';
    else
        error('read_geodata:unknownWCvar',[varLd ' is an '...
            'unknown variable for the WorldClim dataset. Program for this if correct.']);
    end


    if ~regexpbl(wcVar, 'alt')
        fileWcLd = cell(numel(mnthsLd), 1);
        for ii = 1 : numel(mnthsLd)
            if numel(mnthsLd(ii)) == 1 && mnthsLd(ii) < 10
                strMnth = ['0' num2str(mnthsLd(ii))];
            else
                strMnth = num2str(mnthsLd(ii));
            end
            fileTemp = extract_field(dir(fullfile(pathWC, strcat('*', wcVar, '_', strMnth, '.tif' ))), 'name');

            if numel(fileTemp(:)) == 1
                fileWcLd{ii} = fileTemp{1};
            else
                error('readWorldClimTif:wrongFileNumber', [num2str(numel(fileTemp(:))) ...
                    ' files were found. Only one is expected.']);
            end
        end
    else
        fileWcLd = extract_field(dir(fullfile(pathWC, strcat('*', wcVar, '.tif' ))), 'name');
    end

    if isempty(fileWcLd)
        error('read_geodata:wcTifNotFound',['No WorldClim file for month ' ...
            strMnth ' was found.']);
    end
else
    [~, fileWcLd, extLd] = fileparts(pathData);
    fileWcLd = {[fileWcLd, extLd]};
    
%     if regexpbl(fileWcLd, 'prec')
%         wcVar = 'prec';
%     elseif regexpbl(fileWcLd, 'tavg')
%         wcVar = 'tavg';
%     elseif regexpbl(fileWcLd, 'tmin')
%         
%     elseif regexpbl(fileWcLd, 'tmax')
%         
%     elseif regexpbl(fileWcLd, 'alt')
%         
%     else
%         
%     end
end


for ii = 1 : numel(fileWcLd)

    [dataTemp, lonTemp, latTemp, ~] = geoimread(fullfile(pathWC, fileWcLd{ii}));
    
    [data, lon, lat] = crop_grid(dataTemp, lonTemp, latTemp, lonDs, latDs, fr, cropType);
    
    if ii == 1
        sData.(varLd) = nan([numel(fileWcLd), numel(lat), numel(lon)], 'single');
        sData.(varLon) = lon;
        sData.(varLat) = lat;
    end
    sData.(varLd)(ii,:,:) = data;
end
sData.(varLat) = sData.(varLat)(:);

%Set nan values:
sData.(varLd)(sData.(varLd) == -32768) = nan;

%No data values: 32767
valNoData = 32767;
sData.(varLd)(sData.(varLd) == valNoData) = nan;