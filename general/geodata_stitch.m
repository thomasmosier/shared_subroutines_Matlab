function sOut = geodata_stitch(sData1, sData2)

varDate = 'date';
varLon = 'longitude';
varLat = 'latitude';
if iscell(sData1) && iscell(sData2)
    nMod1 = numel(sData1(:));
    nMod2 = numel(sData2(:));
    
    if nMod1 ~= nMod2
       error('geoDataStitch:diffNumMod', ['The first ensemble contains ' ...
           num2tr(nMod1) ' models and the second contains ' ...
           num2str(nMod2) '. The number must be the same.']); 
    end
    
    sOut = cell(nMod1, 1);
    for kk = 1 : nMod1
        if ~isequal(sData1{kk}.(varLon), sData2{kk}.(varLon)) || ~isequal(sData1{kk}.(varLat), sData2{kk}.(varLat))
           warning('geodataStich:gridsDiff', ['The spatial grids for ensemble member ' num2str(kk) ' are not the same. Therefore it is being skipped.']);
           continue
        end
        %Initialize output and combine dates:
        sOut{kk} = sData1{kk};
        sOut{kk}.(varDate) = [sData1{kk}.(varDate); sData2{kk}.(varDate)];
        [sOut{kk}.(varDate), indSrt] = sortrows(sOut{kk}.(varDate));
        blSrt = ~issorted(indSrt);
        
        nTime1 = numel(sData1{kk}.(varDate)(:,1)); 
        nTime2 = numel(sData2{kk}.(varDate)(:,1)); 

        %Find fields to stitch:
        flds1 = fields(sData1{kk});
        flds2 = fields(sData2{kk});

        fldsSame = intersect(flds1, flds2);

        for ii = 1 : numel(fldsSame(:))
            indTime1 = find(nTime1 == size(sData1{kk}.(fldsSame{ii})));
            indTime2 = find(nTime2 == size(sData2{kk}.(fldsSame{ii})));
            if ~iscell(sData1{kk}.(fldsSame{ii})) && ~isempty(indTime1) && ~isempty(indTime2) 
                if indTime1 == indTime2
                    sOut{kk}.(fldsSame{ii}) = cat(indTime1, sData1{kk}.(fldsSame{ii}), sData2{kk}.(fldsSame{ii}));

                    if blSrt
                        nDimsOut = ndims(sOut{kk}.(fldsSame{ii}));
                        if nDimsOut == 3
                            if indTime1 == 1
                                sOut{kk}.(fldsSame{ii}) = sOut{kk}.(fldsSame{ii})(indSrt,:,:);
                            elseif indTime1 == 2
                                sOut{kk}.(fldsSame{ii}) = sOut{kk}.(fldsSame{ii})(:,indSrt,:);
                            else
                                sOut{kk}.(fldsSame{ii}) = sOut{kk}.(fldsSame{ii})(:,:,indSrt);
                            end
                        elseif nDimsOut == 2
                            if indTime1 == 1
                                sOut{kk}.(fldsSame{ii}) = sOut{kk}.(fldsSame{ii})(indSrt,:);
                            else
                                sOut{kk}.(fldsSame{ii}) = sOut{kk}.(fldsSame{ii})(:,indSrt);
                            end
                        end
                    end
                else
                    error('geoDataStitch:diffTimeDim', ['The time dimensions are not the same for field ' fldsSame{ii} '. This has not been programmed for.']); 
                end
            elseif (isempty(indTime1) && ~isempty(indTime2)) || (~isempty(indTime1) && isempty(indTime2))
                if isfield(sOut{kk}, fldsSame{ii})
                   sOut{kk} = rmfield(sOut{kk}, fldsSame{ii}); 
                end
            end
        end
        clear ii
        
        if isfield(sOut{kk}, 'time')
           sOut{kk} = rmfield(sOut{kk}, 'time'); 
        end
    end
    clear kk
% elseif iscell(sData1) && isstruct(sData2)

elseif isstruct(sData1) && isstruct(sData2)
    if ~isequal(sData1.(varLon), sData2.(varLon)) || ~isequal(sData1.(varLat), sData2.(varLat))
       warning('geodataStich:gridsDiff', 'The spatial grids are not the same. Therefore the geodata arrays cannot be stitched.');
       sOut = struct;
       return
    end
    %Initialize output and combine dates:
    sOut = sData1;
    sOut.(varDate) = [sData1.(varDate); sData2.(varDate)];
    [sOut.(varDate), indSrt] = sortrows(sOut.(varDate));
    blSrt = ~issorted(indSrt);

    nTime1 = numel(sData1.(varDate)(:,1)); 
    nTime2 = numel(sData2.(varDate)(:,1)); 

    %Find fields to stitch:
    flds1 = fields(sData1);
    flds2 = fields(sData2);

    fldsSame = intersect(flds1, flds2);

    for ii = 1 : numel(fldsSame(:))
        indTime1 = find(nTime1 == size(sData1.(fldsSame{ii})));
        indTime2 = find(nTime2 == size(sData2.(fldsSame{ii})));
        if ~iscell(sData1.(fldsSame{ii})) && ~isempty(indTime1) && ~isempty(indTime2) 
            if indTime1 == indTime2
                sOut.(fldsSame{ii}) = cat(indTime1, sData1.(fldsSame{ii}), sData2.(fldsSame{ii}));

                if blSrt
                    nDimsOut = ndims(sOut.(fldsSame{ii}));
                    if nDimsOut == 3
                        if indTime1 == 1
                            sOut.(fldsSame{ii}) = sOut.(fldsSame{ii})(indSrt,:,:);
                        elseif indTime1 == 2
                            sOut.(fldsSame{ii}) = sOut.(fldsSame{ii})(:,indSrt,:);
                        else
                            sOut.(fldsSame{ii}) = sOut.(fldsSame{ii})(:,:,indSrt);
                        end
                    elseif nDimsOut == 2
                        if indTime1 == 1
                            sOut.(fldsSame{ii}) = sOut.(fldsSame{ii})(indSrt,:);
                        else
                            sOut.(fldsSame{ii}) = sOut.(fldsSame{ii})(:,indSrt);
                        end
                    end
                end
            else
                error('geoDataStitch:diffTimeDim', ['The time dimensions are not the same for field ' fldsSame{ii} '. This has not been programmed for.']); 
            end
        elseif (isempty(indTime1) && ~isempty(indTime2)) || (~isempty(indTime1) && isempty(indTime2))
            if isfield(sOut, fldsSame{ii})
               sOut = rmfield(sOut, fldsSame{ii}); 
            end
        end
    end
    clear ii

    if isfield(sOut, 'time')
       sOut = rmfield(sOut, 'time'); 
    end
else
    error('geodataOceanNan:unkownGeodataClass', ['The geodata class ' class(sData) ' has not been programmed for.']);
end

