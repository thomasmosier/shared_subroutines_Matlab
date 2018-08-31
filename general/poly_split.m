function [xCell, yCell] = poly_split(x, y)

%See if splits delineated by nans:
indSpltX = find(isnan(x));
indSpltY = find(isnan(y));

if ~isequal(indSpltX, indSpltY)
    error('polySplit:xyMismatch','The splits in the x and y vectors do not match.');
end
crdIn = [x(:), y(:)];

if ~isempty(indSpltX)
    indSpltX = setdiff(indSpltX, [1, numel(x)]);
end

if ~isempty(indSpltX)
    xCell = cell(numel(indSpltX) + 1, 1);
    yCell = xCell;
    
    for ii = 1 : numel(indSpltX) + 1
        if ii == 1
            xCell{1} = crdIn(1:indSpltX(1)-1,1);
            yCell{1} = crdIn(1:indSpltX(1)-1,2);
        elseif ii == numel(indSpltX) + 1
            xCell{end} = crdIn(indSpltX(end)+1:end,1);
            yCell{end} = crdIn(indSpltX(end)+1:end,2);
        else
            xCell{ii} = crdIn(indSpltX(ii-1)+1:indSpltX(ii)-1,1);
            yCell{ii} = crdIn(indSpltX(ii-1)+1:indSpltX(ii)-1,2);
        end
    end
    clear ii
else
    xCell = {crdIn(:, 1)};
    yCell = {crdIn(:, 2)};
end


%Find splits delineated by matching vertices:
for kk = numel(xCell(:)) : -1 : 1
    crdCurr = [xCell{kk}(:), yCell{kk}(:)];
    
    [~,I,~] = unique(crdCurr, 'rows');
    % sort(I)
    ixDup = setdiff(1:numel(crdCurr(:,1)), I);
    crdDup = crdCurr(ixDup,:);

    if numel(ixDup) > 0 %All closed polygons have one duplicate (should be first and last entries)
        %Loop over all duplicate points
        iSplit = []; %This stores the first index of a new segment
        for ii = 1 : numel(ixDup)
            indDup = ismember(crdCurr, crdDup(ii,:), 'rows');
            indDup = find(indDup == 1);
            
            if numel(indDup) ~= 2
                error('polySplit:incorrectDuplicateRows', ['The current vertice is duplicated ' ...
                    num2str(numel(indDup)) ' times. It should have two occurences.']);
            end
            
            if ~isequal(indDup(:), [1; numel(crdCurr(:,1))]) %If duplicates are not first and last, split
                if indDup(1) == 1
                    iSplit = [iSplit, indDup(2)+1];
                else
                    iSplit = [iSplit, indDup(1)];
                end
            end
        end
        
        
        %Add polygon splits to output:
        if ~isempty(iSplit)
            iSplit = sort(iSplit);
        
            nOrg = numel(xCell(:));
            xCell = [xCell; cell(numel(iSplit),1)];
            yCell = [yCell; cell(numel(iSplit),1)];

            xCell{kk} = crdCurr(1:iSplit(1)-1, 1);
            yCell{kk} = crdCurr(1:iSplit(1)-1, 2);
            for ii = 1 : numel(iSplit)
                if ii == numel(iSplit)
                    xCell{nOrg + ii} = crdCurr(iSplit(ii):end, 1);
                    yCell{nOrg + ii} = crdCurr(iSplit(ii):end, 2);
                else
                    xCell{nOrg + ii} = crdCurr(iSplit(ii):iSplit(ii+1)-1, 1);
                    yCell{nOrg + ii} = crdCurr(iSplit(ii):iSplit(ii+1)-1, 2);
                end
            end
        end
    end
end

% %Combine both split delimiters:
% if ~isempty(xCellDup) && ~isempty(xCellNan)
% 	xCell = [xCellNan; xCellDup];
%     yCell = [yCellNan; yCellDup];  
% elseif isempty(xCellDup) && ~isempty(xCellNan)
%     xCell = xCellNan;
%     yCell = yCellNan;    
% elseif ~isempty(xCellDup) && isempty(xCellNan)
%     xCell = xCellDup;
%     yCell = yCellDup;
% else
%     xCell = {x};
%     yCell = {y};
% end
    
    



