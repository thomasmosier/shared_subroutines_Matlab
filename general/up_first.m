function str = up_first(str)

if iscell(str)
    for ii = 1 : numel(str(:))
        str{ii}(1) = upper(str{ii}(1));
    end
else
    str(1) = upper(str(1));
end