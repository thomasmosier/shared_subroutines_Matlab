function sDataIn = struct_2_standard_units(sDataIn, varData, timeStep)

if isfield(sDataIn, ['att' varData])
    strUnits = find_att(sDataIn.(['att' varData]), 'units');
        
    if regexpbl(varData, {'tp','pre','pcp', 'psn','tcw','prcp'}) || strcmpi(varData, 'pr')
        if strcmpi(strUnits, 'm')
            sDataIn.(varData) = sDataIn.(varData)*1000;
            strUnits = 'mm';
        elseif strcmpi(strUnits, 'kg m**-2')
            strUnits = 'mm';
        elseif strcmpi(strUnits, 'kg m-2 s-1')
            scl = 86400;
            if regexpbl(timeStep, {'day', 'daily'})
                sDataIn.(varData) = sDataIn.(varData)*scl;
                strUnits = 'mm/day';
            elseif regexpbl(timeStep, {'mnth', 'month'})
                if numel(size(sDataIn.(varData))) == 3
                    for ii = 1 : numel(sDataIn.date(:,1))
                        sDataIn.(varData)(ii,:,:) = sDataIn.(varData)(ii,:,:)*eomday(sDataIn.date(ii,1), sDataIn.date(ii,2))*scl;
                    end
                    strUnits = 'mm/month';
                elseif numel(size(sDataIn.(varData))) == 2
                    sDataIn.(varData)(:,:) = sDataIn.(varData)(:,:)*eomday(sDataIn.date(1,1), sDataIn.date(1,2))*scl;
                    strUnits = 'mm/month';
                else
                    error('struct_pre_2_mm:unknownSz', ['The data array has size ' num2str(numel(size(sDataIn.(varData)))) ', but only sizes 2 and 3 have been programmed for.']);
                end
                    
                clear ii
            else
                error('struct_pre_2_mm:unknownTimeStep', ['The time step ' timeStep ' has not been programmed for.']);
            end
        elseif strcmpi(strUnits, 'mm h-1')
            if regexpbl(timeStep, {'day', 'daily'})
                sDataIn.(varData) = sDataIn.(varData)*24;
            elseif ~regexpbl(timeStep, {'hr', 'hour'})
                error('struct_pre_2_mm:unknownTimeStep', ['The time step ' timeStep ' has not been programmed for.']);
            end
            strUnits = 'mm';
        elseif ~strcmpi(strUnits, 'mm') && ~((strcmpi(strUnits, 'mm/day') && regexpbl(timeStep, {'day','daily'})) || (strcmpi(strUnits, 'mm/month') && regexpbl(timeStep, 'month')))
            warning('ERA_ds:unknownPreUnits', ['Precipitation units of ' strUnits ...
                ' have not been programmed for.']);
        end
        %Set min to 0:
        sDataIn.(varData)(sDataIn.(varData) < 0) = 0; 
    elseif regexpbl(varData, {'tas','tasmin','tasmax', 'tmp', 'tmn', 'tmx', 'tmean', 'tmin', 'tmax', 't2m', 'mn2t', 'mx2t'}) || strcmpi(varData, 't2')
        if strcmpi(strUnits, 'K') || strcmpi(strUnits, 'Kelvin') 
            sDataIn.(varData) = sDataIn.(varData) - 273.15;
            strUnits = 'Celsius';
        elseif ~(strcmpi(strUnits, 'Celsius') || strcmpi(strUnits, 'C') || strcmpi(strUnits, 'degC') || regexpbl(strUnits, {'deg','celsius'}, 'and'))
            error('ERA_ds:unknownTmpUnits', ['Temperature units of ' strUnits ...
                ' have not been programmed for.']);
        end
    elseif regexpbl(varData, {'rsds', 'rsdt', 'rsdl'})
        if strcmpi(strUnits, 'J m**-2') || strcmpi(strUnits, 'J m-2') 
            sDataIn.(varData) = sDataIn.(varData) / 86400;
            strUnits = 'W m-2';
        elseif ~(strcmpi(strUnits, 'W m-2'))
            error('ERA_ds:unknownRaadUnits', ['Radiation units of ' strUnits ...
                ' have not been programmed for.']);
        end
    elseif regexpbl(varData, {'z','orog'})
        sDataIn.(varData) = squeeze(sDataIn.(varData));

        if regexpbl(strUnits, {'km','kilometer'})
            sDataIn.(varData) = sDataIn.(varData) / 1000;
            strUnits = 'm';
        elseif strcmpi(strUnits, 'm**2 s**-2')
            gp2m = 9.80665; %Constant to transform units of geopotential height into meters
            
            sDataIn.(varData) = sDataIn.(varData) / gp2m; 
            strUnits = 'm';
        elseif ~(strcmpi(strUnits, 'm') || strcmpi(strUnits, 'meters'))
            error('ERA_ds:unknownElevUnits', ['Elevation units of ' strUnits ...
                ' have not been programmed for.']);
        end
    else
        warning('struct_2_standard_units:unknownVar',['The variable ' varData ' has not been programmed for.']);
    end
    
    %Set attributes:
    sDataIn.(['att' varData]) = set_att(sDataIn.(['att' varData]), 'units', strUnits);
else
    sDataIn.(['att' varData]) = {'units', 'unknown'};
    warning('struct_2_standard_units:missingFld',['No field ' char(39) 'att' varData char(39) ' found. This means the units are unknown.'])
end