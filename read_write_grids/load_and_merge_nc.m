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
            nDimCurr = ndims(sDataTemp.(nmsTemp{jj}));

            if nDimCurr == 2 %2D cases
                indTimeCurr = find(size(sDataTemp.(nmsTemp{jj})) == numel(sDataTemp.time));
                if ~isempty(indTimeCurr)
                    sNc.(nmsTemp{jj}) = cat(indTimeCurr, sNc.(nmsTemp{jj}), sDataTemp.(nmsTemp{jj}));
                end               
            elseif nDimCurr == 3 %3D cases (time always first)
                sNc.(nmsTemp{jj}) = cat(1, sNc.(nmsTemp{jj}), sDataTemp.(nmsTemp{jj}));
            end
        end
    end
end

%Turn warning back on (if turned off previously
if blAph == 1
    warning('on', 'ncCal:unknownCal');
end