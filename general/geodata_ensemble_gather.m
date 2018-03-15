function vecOut = geodata_ensemble_gather(sData, var, indLat, indLon)


    
if numel(indLat) ~= numel(indLon)
   error('geodataEnsembleGather:lengthDifferent','The lat and lon indice vectors have different lengths, which is not allowed.');
end

if numel(indLat) ~= 1
    error('geodataEnsembleGather:lengthNot1','The length of the input indice vector is not one. Currently it is mandatory for it to be one. Consider adding functionality.')
end

if iscell(sData)
    if isnumeric(sData{1})
        vecOut = [];
        for ii = 1 : numel(sData(:))
            if ndims(sData{ii}) == 3
                vecOut = [vecOut; squeeze(sData{ii}(:, indLat, indLon))];
            elseif any(size(sData{ii}) == 1)
                vecOut = [vecOut; squeeze(sData{ii}(:))];
            else
                vecOut = [vecOut; squeeze(sData{ii}(indLat, indLon))];
            end
        end
    elseif isstruct(sData{1})
        vecOut = [];
        for ii = 1 : numel(sData(:))
            if ndims(sData{ii}.(var)) == 3
                vecOut = [vecOut; squeeze(sData{ii}.(var)(:, indLat, indLon))];
            elseif any(size(sData{ii}) == 1)
                vecOut = [vecOut; squeeze(sData{ii}.(var)(:))];
            else
                vecOut = [vecOut; squeeze(sData{ii}.(var)(indLat, indLon))];
            end
        end
    end
elseif isstruct(sData)
    if ndims(sData.(var)) == 3
        vecOut = squeeze(sData.(var)(:, indLat, indLon));
    elseif any(size(sData) == 1)
        vecOut = [vecOut; squeeze(sData.(var)(:))];
    else
        vecOut = [vecOut; squeeze(sData.(var)(indLat, indLon))];
    end
else
    error('geodataEnsembleGather:wrongClass',['Input data of class ' class(sData) ' not programmed for.'])
end


