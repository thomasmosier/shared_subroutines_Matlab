function sData = geodata_scale(sData, var, scl)

if iscell(sData)
    for ii = 1 : numel(sData(:))
        sData{ii}.(var) = scl * sData{ii}.(var);
    end
    clear ii
elseif isstruct(sData)
    sData.(var) = scl * sData.(var);
else
    error('geodataScl:unknownClass',['Class ' class ' not recognized as geodata format. Must be cell or struct.'])
end