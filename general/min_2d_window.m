function gridMn = min_2d_window(gridIn, szWind)

if numel(szWind) ~= 1
    error('min_2d_window:szInput','The window must have only one dimension');
end

gridMn = nan(size(gridIn), 'single');

szInput = size(gridIn);

for ii = 1 : szInput(1)
    indRow = [ii-szWind, ii+szWind];
    
    if indRow(1) < 1
        indRow(1) = 1;
    end
    if indRow(2) > szInput(1)
        indRow(2) = szInput(1);
    end
    
    for jj = 1 : szInput(2)
        indCol = [jj-szWind, jj+szWind];
        
        if indCol(1) < 1
            indCol(1) = 1;
        end
        if indCol(2) > szInput(2)
            indCol(2) = szInput(2);
        end
        
        gridMn(ii,jj) = min2d(gridIn(indRow(1):indRow(2), indCol(1):indCol(2)));
    end
end