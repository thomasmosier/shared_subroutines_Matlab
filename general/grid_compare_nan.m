function Zo = grid_compare_nan(Xi,Yi,Zi,Xo,Yo,Zo)

%Method for setting output grid cells to nan:
%Method 1 sets all cells to nan; method 2 only sets exterior cells to
%nan. Both methods take similar times.
indMethod = 1;
    
if any2d(isnan(Zi))
%     figure
%     imagesc(Zi); colorbar; title('Input');
%     figure
%     imagesc(Zo); colorbar; title('Output before nans');

%     tic


    if indMethod == 1
        %Find all nan grid cells in input:
        indNanIn = find(isnan(Zi));
        [rNanIn, cNanIn] = ind2sub(size(Zi), indNanIn);

        xEdgI = box_edg(Xi);
        yEdgI = box_edg(Yi);

        xEdgO = box_edg(Xo);
        yEdgO = box_edg(Yo);

        %Loop over all nan grid cells
        for ll = 1 : numel(indNanIn)
            indLonO = find(xEdgO(2:end)   >= xEdgI(cNanIn(ll))   & xEdgO(1:end-1) <= xEdgI(cNanIn(ll)+1));
            indLatO = find(yEdgO(1:end-1) >= yEdgI(rNanIn(ll)+1) & yEdgO(2:end)   <= yEdgI(rNanIn(ll)));
            
            if ~isempty(indLatO) && ~isempty(indLonO)
                Zo(indLatO(1):indLatO(end), indLonO(1):indLonO(end)) = nan;
            end
        end
        clear ll
    elseif indMethod == 2
        %Looping over longitude slices
        for ii = 1 : numel(Xo(1,:))
            [~, indLoni] = min(abs(Xo(1,ii) - Xi(1,:)));
            indUse = [find(~isnan(Zi(:, indLoni)), 1, 'first'), find(~isnan(Zi(:, indLoni)), 1, 'last')];

            if isempty(indUse)
                Zo(:,ii) = nan;
            else
                Zo(Yo(:,1) > Yi(indUse(1), 1) | Yo(:,1) < Yi(indUse(end), 1), ii) = nan;
            end
        end

        %Loop over latitude slices
        for ii = 1 : numel(Yo(:,1))
            [~, indLati] = min(abs(Yo(ii,1) - Yi(:,1)));
            indUse = [find(~isnan(Zi(indLati,:)), 1, 'first'), find(~isnan(Zi(indLati,:)), 1, 'last')];

            if isempty(indUse)
                Zo(ii,:) = nan;
            else
                Zo(ii, Xo(1,:) < Xi(1,indUse(1)) | Xo(1,:) > Xi(1, indUse(end))) = nan;
            end
        end
    
    end
%     toc
%     figure
%     imagesc(Zo); colorbar; title('Output after nans');
end