function sData = geodata_timestep_convert(sData, var, tStepData, tStepOut, varDate)

varTstep = 'timestep';

if regexpbl(tStepOut, 'month') && regexpbl(tStepData, {'day', 'daily'})
    type = 'avg';
    if strcmpi(var, 'pr') || strcmpi(var, 'pre') 
    %             type = 'sum';
        unitsOut = 'mm/day';
    elseif regexpbl(var, {'tas','dtemp'})
    %             type = 'avg';
        unitsOut = 'Celsius';
    end

    if iscell(var)
        varUnitsIn = var{1};
    else
       varUnitsIn = var; 
    end

    if iscell(sData)
        for kk = 1 : numel(sData(:))
            if regexpbl(tStepData, {'day','daily'})
                sData{kk} = geodata_day2month(sData{kk}, var, type);
                
                sData{kk}.(['att' varUnitsIn]) = set_att(sData{kk}.(['att' varUnitsIn]), 'units', unitsOut);
            elseif regexpbl(tStepData, 'month')
                unitsCurr = find_att(sData{kk}.(['att' varUnitsIn]), 'units');

                if regexpbl(unitsCurr, 'mm/month')
                    nTime = numel(sData{kk}.(varDate));
                    for ll = 1 : nTime
                        sData{kk}.(var)(ll,:,:) = sData{kk}.(var)(ll,:,:) / eomday(sData{kk}.date(ll,1), sData{kk}.date(ll,2));
                    end
                end
            end
            
            sData{kk}.(varTstep) = tStepOut;
        end
        clear kk
    elseif isstruct(sData)
        if regexpbl(tStepData, {'day','daily'})
            sData = geodata_day2month(sData, var, type);
            sData.(varTstep) = tStepOut;

            sData.(['att' varUnitsIn]) = set_att(sData.(['att' varUnitsIn]), 'units', unitsOut);
        elseif regexpbl(tStepData, 'month')
            unitsCurr = find_att(sData.(['att' varUnitsIn]), 'units');

            if regexpbl(unitsCurr, 'mm/month')
                nTime = numel(sData.(varDate));
                for ll = 1 : nTime
                    sData.(var)(ll,:,:) = sData.(var)(ll,:,:) / eomday(sData.date(ll,1), sData.date(ll,2));
                end
            end
        end

        sData.(varTstep) = tStepOut;
    else
        error('geodataTimestepConvert:unknownInputFormat', ['The input format ' class(sData) ' has not been programmed for.']);
    end
elseif ~isequal(tStepOut, tStepData)
    error('geodataTimestepConvert:unknownTimestepCombo', ...
        ['The timesetp combo of input = ' ...
        tStepData ' and desired output = ' tStepOut ...
        ' has not been programmed for.']);
end
    