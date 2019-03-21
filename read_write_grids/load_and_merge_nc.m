function sNc = load_and_merge_nc(filesAll, lonBnds, latBnds, yrBnds, fr, cropType)


if isempty(filesAll(:))
    warning('loadAndMergeNc:noFiles','No file paths were passed to function. This may indicate that no files contained data for requested years.')
    sNc = struct; 
    return
end


%Loop over each of file to load
blAph = 0;
for ii = 1 : length(filesAll(:))
    sDataTemp = read_grid_NC_w_bnds_v2(filesAll{ii}, lonBnds, latBnds, yrBnds, fr, cropType);

    if isempty(sDataTemp.time)
        continue
    end
    
        
    if ii == 1
        sNc = sDataTemp;
    else
        %Convert time to common reference:
        [gcmRefDateTemp, ~] = NC_time_units(sDataTemp.atttime);
        calTemp = NC_cal(sDataTemp.atttime);

        if regexpbl(calTemp, 'unknown')
            calTemp = 'gregorian';
            if ii == 1 
                if ~regexpbl(filesAll{ii}, 'aphro')
                    warning('loadAndMergeNc:unknownCal','The calendar is unknown, but is being set to Gregorian.');
                else 
                    blAph = 1;
                    warning('off', 'ncCal:unknownCal');
                end
            end
        end

        [gcmRefDate, ~] = NC_time_units(sNc.atttime);
        cal =  NC_cal(sDataTemp.atttime);

        if regexpbl(cal, 'unknown')
            cal = 'gregorian';
        end

        dateVecTemp = days_2_date_v2(sDataTemp.time, gcmRefDateTemp, calTemp);
        sDataTemp.time = days_since(gcmRefDate, dateVecTemp, cal);

        %Find all fields to append:
        nmsTemp = fieldnames(sDataTemp);
        nmsNc = fieldnames(sNc);

        if ~isequal(nmsTemp, nmsNc)
            error('loadAndMergeNc:diffNms','The previous and current structure arrays have different fields.');
        end

        for jj = 1 : numel(nmsTemp)
            varCurr = nmsTemp{jj};
            nDimCurr = ndims(sDataTemp.(varCurr));

            if nDimCurr == 2 %2D cases
                indTimeCurr = find(size(sDataTemp.(varCurr)) == numel(sDataTemp.time));
                if ~isempty(indTimeCurr)
                    sNc.(varCurr) = cat(indTimeCurr, sNc.(varCurr), sDataTemp.(varCurr));
                end               
            elseif nDimCurr == 3 %3D cases (time always first)
                szAsn = size(sNc.(varCurr));
                szGet = size(sDataTemp.(varCurr));
                
                if isequal(szAsn(2:3), szGet(2:3)) %Grids are same size
                    sNc.(varCurr) = cat(1, sNc.(varCurr), sDataTemp.(varCurr));
                    
                else %Grids not same size
                    varLd = nmsTemp;
                    %Find latitude variable:
                    varLatTest = {'lat','latitude','row','y','Lat','Latitude','Row','Y','LAT','LATITUDE','ROW'};
                    varLat = blanks(0);
                    for zz = 1 : numel(varLd)
                        if any(strcmpi(varLatTest, varLd{zz}))
                            varLat = varLd{zz};
                            varLd(zz) = []; 
                            break
                        end
                    end

                    if isempty(varLat)
                       error('readGridNcWBnds:noLat',['No latitude variable was found in the NetCDF file being loaded: ' path]) 
                    end

                    %Find longitude variable:
                    varLonTest = {'lon','longitude','col','x','Lon','Longitude','Col','X','LON','LONGITUDE','COL'};
                    varLon = blanks(0);
                    for zz = 1 : numel(varLd)
                        if any(strcmpi(varLonTest, varLd{zz}))
                            varLon = varLd{zz};
                            varLd(zz) = []; 
                            break
                        end
                    end

                    if isempty(varLon)
                       error('readGridNcWBnds:noLon',['No longitude variable was found in the NetCDF file being loaded: ' path]) 
                    end
                    
                    [blDiff, sDataTemp.(varLon), sDataTemp.(varLat), sDataTemp.(varCurr), sNc.(varCurr)] = crd_within(sDataTemp.(varLon), sNc.(varLon), sDataTemp.(varLat), sNc.(varLat), sDataTemp.(varCurr), sNc.(varCurr));
                
                    if blDiff == 0
                        warning('dsLdFlds:grids', ...
                            ['The two grids being stitched are not the same. '...
                            'This will lead to invalid results and probably a script error.']);
                    else
                        sNc.(varLon) = sDataTemp.(varLon); 
                        sNc.(varLat) = sDataTemp.(varLat);
                    end
                    
                    sNc.(varCurr) = cat(1, sNc.(varCurr), sDataTemp.(varCurr));
                end
            end
        end
    end
end

%Turn warning back on (if turned off previously
if blAph == 1
    warning('on', 'ncCal:unknownCal');
end