function [indNm, dimNm] = find_dim(nmsAll, nmFind)


if ~any(strcmp(nmFind,nmsAll))
    indNm = [];
    dimNm = '';
    
    if strcmpi(nmFind, 'latitude')
        varTest = {'lat','latitude','row','y','Lat','Latitude','Row','Y','LAT','LATITUDE','ROW','north_south','south_north'};
    elseif strcmpi(nmFind, 'longitude')
        varTest = {'lon','longitude','col','x','Lon','Longitude','Col','X','LON','LONGITUDE','COL','west_east','east_west'};
    elseif strcmpi(nmFind, 'time')
        varTest = {'time'};
    else
        error('findDim:unknownNm', [nmFind ' is a dimension name that has not been programmed for.'])
    end
    
    for ii = 1 : numel(varTest)
        indNmTemp = find(strcmp(nmsAll, varTest{ii}) == 1);
        if ~isempty(indNmTemp)
           dimNm = varTest{ii};
           indNm = indNmTemp;
           break
        end
    end
    
    if isempty(indNm)
       warning('findDim:noDimFound', [nmFind ' was not matched to any of the available names']); 
    end
else
    indNm = find(strcmp(nmsAll, nmFind) == 1);
    dimNm = nmFind;
end