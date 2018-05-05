function areaOut = area_grid_pt(indLp, lon, lat)

szGrid = [numel(lat), numel(lon)];

areaOut = nan(numel(indLp), 1, 'single');

for ii = 1 : numel(indLp)
    %Convert from indice to row & col:
    [r, c] = ind2sub(szGrid, indLp(ii));

    %Find area of main grid cell that snow is avalanching from:
    if r < szGrid(1) && c < szGrid(2)
        areaOut(ii) = area_geodata(...
            lon(c:c+1)-0.5*diff(lon(c:c+1)), ...
            lat(r:r+1)+0.5*diff(lat(r:r+1)),'e');
    elseif r < szGrid(1) && c == szGrid(2)
        areaOut(ii) = area_geodata(...
            lon(c-1:c)+0.5*diff(lon(c-1:c)), ...
            lat(r:r+1)+0.5*diff(lat(r:r+1)),'e');
    elseif r == szGrid(1) && c < szGrid(2)
        areaOut(ii) = area_geodata(...
            lon(c:c+1)-0.5*diff(lon(c:c+1)), ...
            lat(r-1:r)+1.5*diff(lat(r-1:r)),'e');
    elseif r == szGrid(1) && c == szGrid(2)
        areaOut(ii) = area_geodata(...
            lon(c-1:c)+0.5*diff(lon(c-1:c)), ...
            lat(r-1:r)+1.5*diff(lat(r-1:r)),'e');
    end
end