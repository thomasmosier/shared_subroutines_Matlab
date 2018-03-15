function [attCurr, valNoData] = load_ncatt(pathLd, varCurr)


valNoData = nan;
attCurr = ncinfo(char(pathLd), varCurr);
if ~isempty(attCurr.Attributes)
    attCurr = squeeze(struct2cell(attCurr.Attributes))';

    [rNoVal, ~] = ind2sub(size(attCurr),find(strcmpi(attCurr, 'missing_value') == 1));
    if ~isempty(rNoVal)
        valNoData = attCurr{rNoVal,2};
    end
end
