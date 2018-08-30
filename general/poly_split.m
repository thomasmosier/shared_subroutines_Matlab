function [xCell, yCell] = poly_split(x, y)

%See if splits delineated by nans:
indSpltX = find(isnan(x));
indSpltY = find(isnan(y));

if ~isequal(indSpltX, indSpltY)
    error('polySplit:xyMismatch','The splits in the x and y vectors do not match.');
end

if ~isempty(indSpltX)
    indSpltX = setdiff(indSpltX, [1, numel(x)]);
end

if ~isempty(indSpltX)
    xCellNan = cell(numel(indSpltX) + 1, 1);
    yCellNan = xCellNan;
    
    for ii = 1 : numel(indSpltX) + 1
        if ii == 1
            xCellNan{1} = x(1:indSpltX(1)-1);
            yCellNan{1} = y(1:indSpltX(1)-1);
        elseif ii == numel(indSpltX) + 1
            xCellNan{end} = x(indSpltX(end)+1:end);
            yCellNan{end} = y(indSpltX(end)+1:end);
        else
            xCellNan{ii} = x(indSpltX(ii-1)+1:indSpltX(ii)-1);
            yCellNan{ii} = y(indSpltX(ii-1)+1:indSpltX(ii)-1);
        end
    end
    clear ii
else
    xCellNan = cell(0, 1);
    yCellNan = xCellNan;
end


%Find splits delineated by matching vertices:
crdIn = [x(:), y(:)];
[~,I,~] = unique(crdIn, 'rows');
ixDup = setdiff(1:numel(x), I);

if numel(ixDup) > 0 %All closed polygons have one duplicate (should be first and last entries)
    xCellDup = cell(numel(ixDup), 1);
    yCellDup = xCellDup;
    
    for ii = 1 : numel(ixDup)
        indDup = ismember(crdIn, crdIn(ixDup,:), 'rows');
        indDup = find(indDup == 1);
        if numel(indDup) ~= 2
            error('polySplit:incorrectDuplicateRows', ['The current vertice is duplicated ' ...
                num2str(numel(indDup)) ' times. It should have two occurences.']);
        end
        xCellDup{ii} = crdIn(indDup(1):indDup(2), 1);
        yCellDup{ii} = crdIn(indDup(1):indDup(2), 2);
    end
else
    xCellDup = cell(0, 1);
    yCellDup = xCellDup;
end

%Combine both split delimiters:
if ~isempty(xCellDup) && ~isempty(xCellNan)
	xCell = [xCellNan; xCellDup];
    yCell = [yCellNan; yCellDup];  
elseif isempty(xCellDup) && ~isempty(xCellNan)
    xCell = xCellNan;
    yCell = yCellNan;    
elseif ~isempty(xCellDup) && isempty(xCellNan)
    xCell = xCellDup;
    yCell = yCellDup;
else
    xCell = {x};
    yCell = {y};
end
    
    



