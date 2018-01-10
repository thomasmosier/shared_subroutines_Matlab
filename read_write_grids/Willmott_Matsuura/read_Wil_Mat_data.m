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

function [arrayData, varargout] = read_Wil_Mat_data(dirTS, yrLoad, mnthLoad)   



mnthsUniq = unique(mnthLoad);
nMnths = numel(mnthsUniq);

if numel(mnthLoad) ~= numel(yrLoad) && nMnths == numel(mnthLoad)  
    if numel(yrLoad) == 2
        nYrsLd  = yrLoad(2) - yrLoad(1) + 1;
    else
        nYrsLd = numel(unique(yrLoad));
    end
    
    ldTs = nan(nYrsLd*nMnths,2);
        
    for ii = 1 : nMnths
        ldTs((ii-1)*nYrsLd+ii : ii*nYrsLd,1) = (yrLoad(1):yrLoad(end))';
        ldTs((ii-1)*nYrsLd+ii : ii*nYrsLd,2) = mnthsUniq(ii);
    end
    clear ii
elseif numel(yrLoad) == numel(mnthLoad)
    ldTs = [yrLoad(:), mnthLoad(:)];
else
    error('ASC_file_find:dateFormat','The input data format is unknown.')
end

yrsUniq = unique(ldTs(:,1));

%Find files to load:
filesTS = cell(numel(yrsUniq), 1);
ldTs = nan(numel(yrsUniq), 2);
cntr = 1;
for jj = 1 : numel(yrsUniq)
    %find Willmott and Matsuura file for current year.
    fileCurr = dir([dirTS filesep '*.' num2str(yrsUniq(jj))]);
    
    fileCurr = extract_field(fileCurr, 'name');
    if ~isempty(fileCurr)
        filesTS{cntr} = fileCurr;
        ldTs(jj,:) = [yrsUniq(jj), mnthLoad];
        cntr = cntr + 1;
    end
end
clear jj

if cntr == 1
    warning('readWillMatData:noData', ['No data found for the requested years in folder ' dirTs]);
    return 
end
filesTS(cntr:end) = [];
ldTs(cntr:end,:) = [];


%Initialize 3D numeric array:
arrayData = nan(numel(filesTS(:)), 360 , 720,'single');

cntr = 0;   %Used to assign time indices in output array;

%Willmott and Matsuura data stored by year, so loop over unique years:
fileExist = 0;

for jj = 1 : numel(filesTS(:))
    %Read raw W&M data
    tfid = fopen(fullfile(dirTS,filesTS{jj}));
    tMat = textscan(tfid,'%f');
    tMat = single(reshape(tMat{1}, 14, [])');
    fclose(tfid);

    %%Record subset of months requested for current year:
    %Find months requested for current year:
    mnthsCurr = ldTs(ldTs(:,1) == yrsUniq(jj), 2);
    
    lonWill = tMat(:,1);
    latWill = tMat(:,2);

    indColLon = 1 + (lonWill + 179.75)/0.5;
    indRowLat = 1 + (89.75 - latWill)/0.5;

    for kk = 1 : numel(mnthsCurr)
        for nn = 1 : numel(lonWill)
            arrayData( cntr + kk, indRowLat(nn), indColLon(nn)) = tMat(nn, mnthsCurr(kk) + 2);
        end
    end
       
    %Update counter
    cntr = cntr + numel(mnthsCurr);
end
 
if numel(arrayData(:,1,1)) ~= sum(~isnan(ldTs(:,1)))
    error('read_Wil_Mat_data:nWrong',[num2str(numel(arrayData(:,1,1))) ...
        ' have been loaded from the folder ' char(39) char(dirTS) char(39) ...
        ' but the reqested number to load is ' num2str(numel(ldTs(:,1))) '.']);
end


arrayData(arrayData == -9999) = nan;

if nargout > 1
   varargout{1} = ldTs; 
   
    if nargout > 2
        if fileExist == 1
            varargout{2} = [ length(arrayData(1,1,:)); length(arrayData(1,:,1)); ...
            -180; -90; 0.5; -9999];
        else
           varargout{2} = nan(6,1);
        end
    
        if nargout > 3
            if fileExist == 1
                varargout{3} = {'NCOLS';'NROWS';'XLLCORNER';'YLLCORNER';'CELLSIZE';'NODATA_VALUE'};
            else
                varargout{3} = cell(6,1);
            end
        end
   end
end
