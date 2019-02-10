function sData = load_NC_ensemble(pathLd, varLd, varOut, lonUse, latUse, yrsLd, fr, crop, varargin)

blWgt = 0;
lonWgt = [];
latWgt = [];
timeStepUse = [];
varTStep = [];
timeCombine = [];
for ii = 1 : numel(varargin(:))
    if regexpbl(varargin{ii}, {'area', 'wgt'}, 'and')
        blWgt = 1;
        lonWgt = varargin{ii+1};
        latWgt = varargin{ii+2};
    elseif regexpbl(varargin{ii}, {'time', 'convert'}, 'and')
        timeStepUse = varargin{ii+1};
        varTStep  = varargin{ii+2};
    elseif regexpbl(varargin{ii}, {'time', 'combine'}, 'and')
        timeCombine = 1;
    end
end
clear ii

varDate = 'date';

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

%Loop over input paths
for kk = 1 : nMem
    if iscell(pathLd{1}) %Check first ensemble member here
        sDataTemp = cell(numel(varLd),1);

        for ll = 1 : numel(varLd(:))
            disp(['Loading ensemble member ' num2str(kk) ' of ' num2str(nMem) ' (variable ' num2str(ll) ' of ' num2str(numel(varLd(:))) ').']);

            sDataTemp{ll} = read_geodata_v2(pathLd{ll}{kk}, varLd{ll}, lonUse, latUse, yrsLd, fr, crop, 'units', 'standard', 'no_disp');
            sDataTemp{ll} = struct_2_standard_units(sDataTemp{ll}, varLd{ll}, sDataTemp{ll}.timestep);
            
            if ~isempty(timeStepUse)
                sDataTemp{ll} = geodata_timestep_convert(sDataTemp{ll}, varLd{ll}, sDataTemp{ll}.(varTStep), timeStepUse, varDate);
            end
        end
        clear ll

        if regexpbl(varOut, 'dtemp')
            sData{kk} = sDataTemp{1};
            if regexpbl(varLd{1}, 'tasmax') && regexpbl(varLd{2}, 'tasmin') 
                sData{kk}.(varOut) = sDataTemp{1}.(varLd{1}) - sDataTemp{2}.(varLd{2});
            elseif regexpbl(varLd{2}, 'tasmax') && regexpbl(varLd{1}, 'tasmin') 
                sData{kk}.(varOut) = sDataTemp{2}.(varLd{2}) - sDataTemp{1}.(varLd{1});
            else
                error('loadNcEnsemble:unknownDtempVar','The variables programmed for dtemp do not seem correct.')
            end
        else
            error('loadNcEnsemble:unknownVar',['Variable ' varOut ' is not expected to require loading multiple variables.']);
        end
    else
        if iscell(varLd)
            if numel(varLd(:)) == 1
                varLd = varLd{1};
            else
                error('loadNcEnsemble:multVar','Multiple input variables selected to load, but format of file paths is not correct to support loading multiple variables.');
            end
        end
        
        disp(['Loading ensemble member ' num2str(kk) ' of ' num2str(nMem) '.']);

        sData{kk} = read_geodata_v2(pathLd{kk}, varLd, lonUse, latUse, yrsLd, fr, crop, 'units', 'standard', 'no_disp');
        sData{kk} = struct_2_standard_units(sData{kk}, varLd, sData{kk}.timestep);
        
        if ~isequal(varLd, varOut)
            sData{kk}.(varOut) = sData{kk}.(varLd);
            sData{kk} = rmfield(sData{kk}, varLd);
        end
        
        if ~isempty(timeStepUse)
            sData{kk} = geodata_timestep_convert(sData{kk}, varOut, sData{kk}.(varTStep), timeStepUse, varDate);
        end
    end
    
    if blWgt == 1
        varLon = 'longitude';
        varLat = 'latitude';

        sData{kk}.(varOut) = area_conserve_remap(sData{kk}.(varLon), sData{kk}.(varLat), sData{kk}.(varOut), lonWgt, latWgt);

        sData{kk}.(varLon) = lonWgt;
        sData{kk}.(varLat) = latWgt;
    end
end
clear kk


%If requested, combine input data along time dimension
varTime = 'time';
varDate = 'date';
varTimeBnds = 'time_bnds';

if timeCombine
    sTemp = sData{1};
    
    if numel(sData(:)) > 1
        for ii = 2 : numel(sData(:))
            sTemp.(varOut) = cat(1, sTemp.(varOut), sData.(varOut));
            sTemp.(varTime) = cat(1, sTemp.(varTime), sData.(varTime));
            sTemp.(varDate) = cat(1, sTemp.(varDate), sData.(varDate));
            if isfield(sTemp, varTimeBnds)
                sTemp.(varTimeBnds) = cat(1, sTemp.(varTimeBnds), sData.(varTimeBnds));
            end
        end
        
        sData = sTemp;
        clear sTemp
    end

end