function sData = load_NC_ensemble(pathLd, varLd, varOut, lonUse, latUse, yrsLd, fr, crop)

if ischar(pathLd)
    pathLd = {pathLd};
end

%Determine number of ensemble members for each simulation:
if iscell(pathLd{1})
    nMem = numel(pathLd{1});
else
    nMem = numel(pathLd);
end

%Cell/structure arrays for main data to be processed:
sData = cell(nMem, 1);

for kk = 1 : nMem
    if iscell(pathLd{1})
        sDataTemp = cell(numel(varLd));

        for ll = 1 : numel(varLd(:))
            sDataTemp{ll} = read_geodata_v2(pathLd{ll}{kk}, varLd{ll}, lonUse, latUse, yrsLd, fr, crop, 'units', 'standard');
            sDataTemp{ll} = struct_2_standard_units(sDataTemp{ll}, varLd{ll}, sDataTemp{ll}.timestep);
        end
        clear ll

        if regexpbl(varOut, 'dtemp')
            sData{kk} = sDataTemp{1};
            if regexpbl(varLd{1}, 'tasmax') && regexpbl(varLd{2}, 'tasmin') 
                sData{kk}.(varOut) = sDataTemp{1}.(varLd{1}) - sDataTemp{2}.(varLd{2});
            elseif regexpbl(varLd{2}, 'tasmax') && regexpbl(varLd{1}, 'tasmin') 
                sData{kk}.(varOut) = sDataTemp{2}.(varLd{2}) - sDataTemp{1}.(varLd{1});
            else
                error('GISSCompare:unknownDtempVar','The variables programmed for dtemp do not seem correct.')
            end
        else
            error('GISSCompare:unknownVar',['Variable ' varOut ' is not expected to require loading multiple variables.']);
        end
    else
        sData{kk} = read_geodata_v2(pathLd{kk}, varOut, lonUse, latUse, yrsLd, fr, crop, 'units', 'standard');
        sData{kk} = struct_2_standard_units(sData{kk}, varOut, sData{kk}.timestep);
    end
end
clear kk