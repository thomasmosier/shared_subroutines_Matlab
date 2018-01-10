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
function [foldMainOut, foldStdOut] = season(foldTs,sMeta,calcTyp,varargin)
%This function is converts monthly time-series to seasonal time-series.

%%INPUTS:
%sMeta.mnths = List of months to use in calculating season
%'strType': Either 'sum' (sums data within month vector) or 'mean' (averages data within month vector).
%varargin: first argument used to set scale parameter (default is 1)

if ~isempty(varargin) && ~isempty(varargin{1})
    sclTs = varargin{1};
else
    sclTs = 1;
end

%In future version, change so that it can produce seasonal from ascii or netCDF
extTyp = 'asc';

%Find fiels in folder specified:
gridFiles = find_files(foldTs,extTyp);

%%Create Output directory:
dirOutput = mk_output_dir(fullfile(foldTs, gridFiles{1}), 'seasonal');
disp(['Seasonal output processing has begun.' char(10) ...
    'The seasonal output files will be written to: ' dirOutput]);

%Find indices seperating year and month (should be in format '*_year_mth.asc')
undInd = regexpi(gridFiles{1},'_');
if length(undInd) > 2
   undInd = undInd(end-1:end);
elseif length(undInd) < 2
    error('climFromTs:und',['File format does not contain expected number '...
        'of underscores.']);
end

if sclTs ~= 1
    disp(['Warning: A scale factor of ' num2str(sclTs) ' is being used.'])
end

%%Find all sets of years present in the folder:
unqYrs = cell(0,1);
for ii = 1 : length(gridFiles(:))
    currYr = gridFiles{ii}(undInd(1)+1:undInd(2)-1);
    if sum(strcmpi(unqYrs, char(currYr))) == 0
        unqYrs{end+1} = currYr;
    end
end

%Loop over all unique years in folder:
for kk = 1 : length(unqYrs(:))    
    %Initialize list of files to use:
    filesUseBl = zeros(length(gridFiles(:)),1);

    %Find all files in directory that match current year and each of the
    %months:
    for ii = 1 : length(sMeta.mnths(:)) 
        filesCurr = cellfun(@(x) ...
            strcmpi(x(undInd(1)+1:undInd(2)-1),unqYrs{kk}) ...
            && str2double(x(undInd(2) + 1 : end - 4)) == sMeta.mnths(ii) ...
            , gridFiles);

        filesUseBl(filesCurr) = 1;
    end

    %Check that correct number of files selected:
    if sum(filesUseBl) ~= length(sMeta.mnths(:))
       error('climFromTs:incomYrs',['An incorrect number of time ' ...
           'series files were selected from ' foldTs '.' ]); 
    end

    %Select files to use in current iteration:
    filesUse = gridFiles(filesUseBl == 1);

    %Load files and calculate seasonal average:
    for jj = 1 : length(filesUse(:))
        [dataTsCurr, hdrCurr] = read_ESRI(fullfile(foldTs, filesUse{jj}));

        if jj == 1
            dataTs = zeros([length(sMeta.mnths(:)),size(dataTsCurr)],'single');
            hdrFir = hdrCurr;
        else
            if ~isequal(hdrFir, hdrCurr)
                error('climFromTs:misalignHdr',['The header ' ...
                    'corresponding to ' filesUse{jj} ' does not align ' ...
                    'with that of ' filesUse{1} '.']);
            end
        end

        dataTs(jj,:,:) = dataTsCurr * sclTs;
    end


    %Create output file name:
    nmAvg = char(filesUse{1}(1:undInd(2)));
    for ll = 1 : length(sMeta.mnths(:))
        nmAvg = [nmAvg num2str(sMeta.mnths(ll)) '.'];
    end
    nmAvg = [nmAvg 'asc'];

    %Write sum or average to file:
    if regexpbl(calcTyp,'sum')
        dataSum = single(squeeze(nansum(dataTs,1)));
        %Write seasonal sum to file:
        %Determine how many decimals to use:
        decSum = dec_write(dataSum);
        nmSum = [nmAvg(1:end-4), '_sum', nmAvg(end-3:end)];
        foldMainOut = fullfile(dirOutput, 'sum');
        fullPathSum = fullfile(foldMainOut, nmSum);
        write_ESRI_v4(dataSum, hdrFir, fullPathSum, decSum);
    elseif regexpbl(calcTyp,{'mean','avg','average'})
        dataAvg = squeeze(nanmean(dataTs,1));
        %Determine how many decimals to use:
        decAvg = dec_write(dataAvg);
        foldMainOut = fullfile(dirOutput, 'avg');
        fullPathClim = fullfile(foldMainOut, [nmAvg(1:end-4), '_avg', nmAvg(end-3:end)]);
        write_ESRI_v4(dataAvg, hdrFir, fullPathClim, decAvg);
    end

    %Write standard deviation:
    dataStd = single(squeeze(nanstd(dataTs,1)));
    %Determine how many decimals to use:
    decStd = dec_write(dataStd);
    %Path:
    nmStd = [nmAvg(1:end-4), '_std', nmAvg(end-3:end)];
    foldStdOut = fullfile(dirOutput, 'std');
    fullPathStd = fullfile(foldStdOut, nmStd);
    write_ESRI_v4(dataStd, hdrFir, fullPathStd, decStd);

    %Display update:
    fracComplt = round2(100*kk/length(unqYrs(:)),2);
    warning('off','MATLAB:gui:latexsup:UnableToInterpretTeXString');
    display(['The seasonal average ' nmAvg ...
        ' has been written.' char(10) 'The season function is '...
        num2str(fracComplt) '% complete']);
    warning('on','MATLAB:gui:latexsup:UnableToInterpretTeXString');
end
