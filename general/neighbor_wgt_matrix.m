function wgtMat = neighbor_wgt_matrix(szMat, strType)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nMat = szMat(1)*szMat(2);

%Initialize
% wgtMat = sparse(nMat, nMat);

indCor = unique(sub2ind(szMat, [1,1, szMat(1), szMat(1)], [1, szMat(2), 1, szMat(2)]));

indEdg = unique(setdiff(sub2ind(szMat, ...
    [ones(1,szMat(2)), szMat(1)*ones(1,szMat(2)), (1:szMat(1)), (1:szMat(1))], ...
    [(1:szMat(2)), (1:szMat(2)), ones(1,szMat(1)), szMat(2)*ones(1,szMat(1))]), ...
    indCor));

indCnt = setdiff((1:szMat(1)*szMat(2)), unique([indCor(:);indEdg(:)]));

if regexpbl(strType, 'all')
    [ic, icd] = ixneighbors(ones(szMat,'single'));
    
    %Weights for corners, edges, and interior grid cells
    val = [1/3, 1/5, 1/8];
    
%     wgtMat(sub2ind([nMat, nMat], ic(indCorO), icd(indCorO))) = 1 / 3;
%     wgtMat(sub2ind([nMat, nMat], ic(indEdgO), icd(indEdgO))) = 1 / 5;
%     wgtMat(sub2ind([nMat, nMat], ic(indCntO), icd(indCntO))) = 1 / 8;

elseif regexpbl(strType, 'perp')
    [ic, icd] = ixneighbors(ones(szMat,'single'), '', 4);
%     
%     indCorO = ismember(ic, indCor);
%     indEdgO = ismember(ic, indEdg);
%     indCntO = ismember(ic, indCnt);
    
    %Weights for corners, edges, and interior grid cells
    val = [1/2, 1/3, 1/4];
    
%     wgtMat = sparse(...
%         [ic(indCorO); ic(indEdgO); ic(indCntO)], ...
%         [icd(indCorO); icd(indEdgO); icd(indCntO)], ...
%         [(1/2)*ones(sum(indCorO), 1); (1/3)*ones(sum(indEdgO), 1); (1/4)*ones(sum(indCntO), 1)], ...
%         nMat, nMat);
end

indCorO = ismember(ic, indCor);
indEdgO = ismember(ic, indEdg);
indCntO = ismember(ic, indCnt);

wgtMat = sparse(...
    [ic(indCorO); ic(indEdgO); ic(indCntO)], ...
    [icd(indCorO); icd(indEdgO); icd(indCntO)], ...
    [val(1)*ones(sum(indCorO), 1); val(2)*ones(sum(indEdgO), 1); val(3)*ones(sum(indCntO), 1)], ...
    nMat, nMat);

% wgtMat2 = sparse(nMat, nMat);
% for ii = 1 : szMat(1)
%     for jj = 1 : szMat(2)
%         if regexpbl(strType, 'all')
%             if ii == 1 && jj == 1
%                 neighbors = sub2ind(szMat, ii+[0, 0, 1, 1], jj+[0, 1, 0, 1]);
%             elseif ii == 1 && jj == szMat(2)
%                 neighbors = sub2ind(szMat, ii+[0, 0, 1, 1], jj+[0, -1, 0, -1]);
%             elseif ii == szMat(1) && jj == 1
%                 neighbors = sub2ind(szMat, ii+[0, 0, -1, -1], jj+[0, 1, 0, 1]);
%             elseif ii == szMat(1) && jj == szMat(2)
%                 neighbors = sub2ind(szMat, ii+[0, 0, -1, -1], jj+[0, -1, 0, -1]);
%             elseif ii == 1
%                 neighbors = sub2ind(szMat, ii+[0, 0, 0, 1, 1, 1], jj+[-1, 0, 1, -1, 0, 1]);
%             elseif jj == 1
%                 neighbors = sub2ind(szMat, ii+[-1, 0, 1, -1, 0, 1], jj+[0, 0, 0, 1, 1, 1]);
%             elseif ii == szMat(1) && jj == szMat(2)
%                 neighbors = sub2ind(szMat, ii+[0, 0, -1, -1], jj+[0, -1, 0, -1]);
%             elseif ii == szMat(1) 
%                 neighbors = sub2ind(szMat, ii+[0, 0, 0, -1, -1, -1], jj+[-1, 0, 1, -1, 0, 1]);
%             elseif jj == szMat(2)   
%                 neighbors = sub2ind(szMat, ii+[-1, 0, 1, -1, 0, 1], jj+[0, 0, 0, -1, -1, -1]);
%             else
%                 neighbors = sub2ind(szMat, ii+[-1, -1, -1, 0, 0, 0, 1, 1, 1], jj+[-1, 0, 1, -1, 0, 1, -1, 0, 1]);
%             end
%         elseif regexpbl(strType, 'perp')
%             if ii == 1 && jj == 1
%                 neighbors = sub2ind(szMat, ii+[0, 1], jj+[1, 0]);
%             elseif ii == 1 && jj == szMat(2)
%                 neighbors = sub2ind(szMat, ii+[0, 1], jj+[-1, 0]);
%             elseif ii == szMat(1) && jj == 1
%                 neighbors = sub2ind(szMat, ii+[0, -1], jj+[1, 0]);
%             elseif ii == szMat(1) && jj == szMat(2)
%                 neighbors = sub2ind(szMat, ii+[0, -1], jj+[-1, 0]);
%             elseif ii == 1
%                 neighbors = sub2ind(szMat, ii+[0, 1, 0], jj+[-1, 0, 1]);
%             elseif jj == 1
%                 neighbors = sub2ind(szMat, ii+[-1, 0, 1], jj+[0, 1, 0]);
%             elseif ii == szMat(1) 
%                 neighbors = sub2ind(szMat, ii+[0, -1, 0], jj+[-1, 0, 1]);
%             elseif jj == szMat(2)   
%                 neighbors = sub2ind(szMat, ii+[-1, 0, 1], jj+[0, -1, 0]);
%             else
%                 neighbors = sub2ind(szMat, ii+[-1, -1, -1, 0, 0, 0, 1, 1, 1], jj+[-1, 0, 1, -1, 0, 1, -1, 0, 1]);
%             end
%         end
%         
%         target = sub2ind(szMat, ii, jj);
%         wgtMat2(target, neighbors) = 1 / numel(neighbors);
%     end
% end
