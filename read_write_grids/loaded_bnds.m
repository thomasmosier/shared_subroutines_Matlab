function [lonBnds, latBnds] = loaded_bnds(lon, lat)


errThresh = 0.01;

if any(size(lon) == 1) && any(size(lat) == 1)
    stepLon = mean(abs(diff(lon)));
    stepLat = mean(abs(diff(lat)));
    
    if lat(1) > lat(end)
        latBnds = [lat(1) + (0.5+errThresh)*stepLat, lat(end) - (0.5+errThresh)*stepLat];
    else
        latBnds = [lat(1) - (0.5+errThresh)*stepLat, lat(end) + (0.5+errThresh)*stepLat];
    end

    if lon(1) < lon(end)
        lonBnds = [lon(1) - (0.5+errThresh)*stepLon, lon(end) + (0.5+errThresh)*stepLon];
    else
        lonBnds = [lon(1) + (0.5+errThresh)*stepLon, lon(end) - (0.5+errThresh)*stepLon];
    end
else
    stepLonTemp = numel(lon(:,1));
    for ii = 1 : numel(lon(:,1))
        stepLonTemp(ii) = mean(abs(diff(lon(ii,:))));
    end
    stepLon = nanmean(stepLonTemp);
   
    stepLatTemp = numel(lat(1,:));
    for ii = 1 : numel(lat(1,:))
        stepLatTemp(ii) = mean(abs(diff(lat(:,ii))));
    end
    stepLat = nanmean(stepLatTemp);
    
    
    if lat(1,1) > lat(end,1)
        latBnds = [max2d(lat) + (0.5+errThresh)*stepLat, min2d(lat) - (0.5+errThresh)*stepLat];
    else
        latBnds = [min2d(lat) - (0.5+errThresh)*stepLat, max2d(lat) + (0.5+errThresh)*stepLat];
    end

    if lon(1,1) < lon(1,end)
        lonBnds = [min2d(lon) - (0.5+errThresh)*stepLon, max2d(lon) + (0.5+errThresh)*stepLon];
    else
        lonBnds = [max2d(lon) + (0.5+errThresh)*stepLon, min2d(lon) - (0.5+errThresh)*stepLon];
    end
end

lonBnds = sort(lonBnds);
latBnds = sort(latBnds);

