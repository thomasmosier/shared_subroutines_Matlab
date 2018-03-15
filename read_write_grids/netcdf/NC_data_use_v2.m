% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function [gcmTs, attData, attGen, nmData] = NC_data_use_v2(path, nmData, tInd, latInd, lonInd)



%Get dimension names from netCDf file:
dimRaw = ncinfo(path,nmData);
dimNC = cell(numel(dimRaw.Dimensions(:)), 1);
for ii = 1 : numel(dimRaw.Dimensions(:))
    dimNC{ii} = dimRaw.Dimensions(ii).Name;
end

%Define desired order:
indD = {tInd, latInd, lonInd};
dimOrdered = {'time';'latitude';'longitude'};
    
if ~isempty(dimNC) 
    if ~any(strcmp(dimOrdered{2},dimNC))
        dimOrdered{2} = 'Latitude';  
    end
    if ~any(strcmp(dimOrdered{2},dimNC))
        dimOrdered{2} = 'lat';  
    end
    if ~any(strcmp(dimOrdered{2},dimNC))
       error('NC_data_use:latNm','Cannot find correct name of latitude dimension.');
    end
    if ~any(strcmp(dimOrdered{3},dimNC))
        dimOrdered{3} = 'Longitude';  
    end
    if ~any(strcmp(dimOrdered{3},dimNC))
        dimOrdered{3} = 'lon';  
    end
    if ~any(strcmp(dimOrdered{3},dimNC))
       error('NC_data_use:lonNm','Cannot find correct name of longitude dimension.');
    end

    sortDim = (1:3);
    %If dimensions not ordered, reorder:
    if ~all(strcmp(dimOrdered,dimNC))
        for ii = 1 : numel(dimOrdered)
            for jj = 1 : numel(dimNC)
                if strcmpi(dimOrdered{ii},dimNC{jj})
                    sortDim(ii) = jj;
                end
            end
        end
    end
else %isempty(dimNC)
    error('NC_data_use:wrongVar','The variable provided the function is not correct.');
end


dataSz = dimRaw.Size;
blDimOrd = 1;
for ii = 1 : numel(dataSz);
    if dataSz(ii) < indD{sortDim(ii)}(end)
       blDimOrd = 0; 
    end
end


if blDimOrd == 0 %Dimensions must be out of order! Caused by person who made NetCDF file... :(
	ncid = netcdf.open(path, 'NC_NOWRITE');
    blNc = 0;
    cntrNc = 0;
    varNc = cell(0,1);
    while blNc == 0
        try
            [varname, ~, ~, ~] = netcdf.inqVar(ncid, cntrNc);
            varNc{cntrNc+1} = varname;
        catch %'read_geodata:splice'
            blNc = 1;
        end
        cntrNc = cntrNc + 1;
    end
    netcdf.close(ncid);
    
    dimNC = varNc;
    
    if ~any(strcmp(dimOrdered{2},dimNC))
        dimOrdered{2} = 'Latitude';  
    end
    if ~any(strcmp(dimOrdered{2},dimNC))
        dimOrdered{2} = 'lat';  
    end
    if ~any(strcmp(dimOrdered{2},dimNC))
       error('NC_data_use:latNm','Cannot find correct name of latitude dimension.');
    end
    if ~any(strcmp(dimOrdered{3},dimNC))
        dimOrdered{3} = 'Longitude';  
    end
    if ~any(strcmp(dimOrdered{3},dimNC))
        dimOrdered{3} = 'lon';  
    end
    if ~any(strcmp(dimOrdered{3},dimNC))
       error('NC_data_use:lonNm','Cannot find correct name of longitude dimension.');
    end
    
    dataSz = ncinfo(path, nmData);
    dataSz = dataSz.Size;
    
    varSz = nan(3,1);
    for ii = 1 : 3
        dataSzTemp = ncinfo(path, dimOrdered{ii});

        varSz(ii) = dataSzTemp.Size;
    end
    
    
    %If dimensions not ordered, reorder:
    if ~isequal(varSz, dataSz)
        for ii = 1 : numel(dataSz)
            for jj = 1 : numel(varSz)
                if dataSz(ii) == varSz(jj)
                    sortDim(ii) = jj;
                end
            end
        end
    end
end



%If data out of order, must load in pieces:
ld = zeros(3,1);
if ~issorted(indD{sortDim(1)}) || any(diff(indD{sortDim(1)}) ~= 1)
    ld(1) = 1;
end
if ~issorted(indD{sortDim(2)}) || any(diff(indD{sortDim(2)}) ~= 1)
    ld(2) = 1;
end
if ~issorted(indD{sortDim(3)}) || any(diff(indD{sortDim(3)}) ~= 1)
    ld(3) = 1;
end


if sum(ld) == 0 %If no pieces need to be re-ordered:
    gcmTs = single(ncread(path, nmData, [indD{sortDim(1)}(1), indD{sortDim(2)}(1), indD{sortDim(3)}(1)], [numel(indD{sortDim(1)}), numel(indD{sortDim(2)}), numel(indD{sortDim(3)})]));
elseif sum(ld) == 1
    gcmTs = zeros([numel(indD{sortDim(1)}), numel(indD{sortDim(2)}), numel(indD{sortDim(3)})],'single');
    if ld(1) == 1
        [~, indSort1] = sort(indD{sortDim(1)});
        if all(diff(indSort1) ==1)
            try
                gcmTs = single(ncread(path, nmData, [indD{sortDim(1)}(1), indD{sortDim(2)}(1), indD{sortDim(3)}(1)], [numel(indD{sortDim(1)}), numel(indD{sortDim(2)}), numel(indD{sortDim(3)})]));
                catch ME
                    if (strcmp(ME.identifier,'MATLAB:nomem'))
                        indX = round(0.5*indD{sortDim(1)}(1));
                        gcmTs(1:indX,:,:) = single(ncread(path, nmData, [indD{sortDim(1)}(1), indD{sortDim(2)}(1), indD{sortDim(3)}(1)], [indX, numel(indD{sortDim(2)}), numel(indD{sortDim(3)})]));
                        gcmTs(indX:end,:,:) = single(ncread(path, nmData, [indD{sortDim(1)}(1), indD{sortDim(2)}(1), indD{sortDim(3)}(1)], [numel(indD{sortDim(1)}) - indX, numel(indD{sortDim(2)}), numel(indD{sortDim(3)})]));
                    end
            end
            gcmTs = gcmTs(indSort1,:,:);
        else
            for ii = 1 : numel(indD{sortDim(1)})
                gcmTs(ii,:,:) = single(ncread(path, nmData, [indD{sortDim(1)}(ii), indD{sortDim(2)}(1), indD{sortDim(3)}(1)], [1, numel(indD{sortDim(2)}), numel(indD{sortDim(3)})]));
            end
        end
    elseif ld(2) == 1
        for ii = 1 : numel(indD{sortDim(2)})
            gcmTs(:,ii,:) = single(ncread(path, nmData, [indD{sortDim(1)}(1), indD{sortDim(2)}(ii), indD{sortDim(3)}(1)], [numel(indD{sortDim(1)}), 1, numel(indD{sortDim(3)})]));
        end
    elseif ld(3) == 1
        for ii = 1 : numel(indD{sortDim(3)})
            gcmTs(:,:,ii) = single(ncread(path, nmData, [indD{sortDim(1)}(1), indD{sortDim(2)}(1), indD{sortDim(3)}(ii)], [numel(indD{sortDim(1)}), numel(indD{sortDim(2)}), 1]));
        end
    end
    
elseif sum(ld) > 1
    gcmTs = zeros([numel(indD{1}), numel(indD{2}), numel(indD{3})],'single');
    warning('NC_data_use:d2sep',[num2str(sum(ld)) ' dimensions are ' ...
        'out of order. This increases NetCDF read time.']);
    for ii = 1 : numel(indD{1})
        for jj = 1 : numel(indD{2})
            for kk = 1 : numel(indD{3})
                vecCurr = [indD{1}(ii), indD{2}(jj), indD{3}(kk)];
                gcmTs(ii,jj,kk) = single(ncread(path, nmData, vecCurr(sortDim), [1, 1, 1]));
            end
        end
    end
    sortDim = (1:3);
end
    

%Reorder dimenions if neccesarry:
if ~issorted(sortDim)
    [~, indSort] = sort(sortDim);
    gcmTs = permute(gcmTs, indSort);
end

%Load general attributes:
attGen = ncinfo(path);
attGen = squeeze(struct2cell(attGen.Attributes))';

%Load data attributes:
attData = ncinfo(path,nmData);
attData = squeeze(struct2cell(attData.Attributes))';

end