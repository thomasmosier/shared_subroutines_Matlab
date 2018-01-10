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

function var_consistency_chk(sMeta,sPath)



for ii = 1 : length(sMeta.metVar)
    if ~isempty(regexpi(sMeta.metVar{ii},'pre'))
        varList = {sMeta.metVar{ii},'pr'};
    elseif ~isempty(regexpi(sMeta.metVar{ii},'tmn'))
        varList = {sMeta.metVar{ii},'tmin','tasmin'};
    elseif ~isempty(regexpi(sMeta.metVar{ii},'tmx'))
        varList = {sMeta.metVar{ii},'tmax','tasmax'};
    elseif ~isempty(regexpi(sMeta.metVar{ii},'tmp'))
        varList = {sMeta.metVar{ii},'tmean','tas'};
    else
        error('var_consistency_chk:unknownVar',[sMeta.metVar{ii} ...
            ' is an unknown variable name.']);
    end
    
    %Check historical time-series
    if ~isempty(sPath.hisObs{ii})
        varhisTS = 0;
        for jj = 1 : length(varList)
            if ~isempty(regexpi(sPath.hisObs{ii},varList{jj}));
                varhisTS = 1;
            end
        end

        if varhisTS == 0
            warning('var_consistency_chk:diffVar',['The data file or '...
                'directory chosen is ' char(39) sPath.hisObs{ii} char(39) ...
                ' for the variable ' char(39) sMeta.metVar{ii} char(39) ...
                char(10) '.  The automated check believes these may be inconsistent.']);
        end
    end
    
    %Check historical GCM
    if ~isempty(sPath.hisMod{ii})
        varhisTS = 0;
        for jj = 1 : length(varList)
            if ~isempty(regexpi(sPath.hisMod{ii},varList{jj}));
                varhisTS = 1;
            elseif regexpbl(sPath.hisMod{ii}, 'gpcc') && regexpbl(varList{jj}, 'pre')
                varhisTS = 1;
            end
        end

        if varhisTS == 0
            warning('var_consistency_chk:diffVar',['The data file or '...
                'directory chosen is ' char(39) sPath.hisMod{ii} char(39) ...
                ' for the variable ' char(39) sMeta.metVar{ii} char(39) ...
                char(10) '.  The automated check believes these may be inconsistent.']);
        end
    end
    
    %Check high-resolution climatology
    if ~isempty(sPath.hrClim{ii})
        varhisTS = 0;
        for jj = 1 : length(varList)
            if ~isempty(regexpi(sPath.hrClim{ii},varList{jj}));
                varhisTS = 1;
            end
        end

        if varhisTS == 0
            warning('var_consistency_chk:diffVar',['The data file or '...
                'directory chosen is ' char(39) sPath.hrClim{ii} char(39) ...
                ' for the variable ' char(39) sMeta.metVar{ii} char(39) ...
                char(10) '.  The automated check believes these may be inconsistent.']);
        end
    end
    
    %Check projected Time-series
    if ~isempty(sPath.projMod{ii})
        varhisTS = 0;
        for jj = 1 : length(varList)
            if ~isempty(regexpi(sPath.projMod{ii},varList{jj}));
                varhisTS = 1;
            end
        end

        if varhisTS == 0
            warning('var_consistency_chk:diffVar',['The data file or '...
                'directory chosen is ' char(39) sPath.projMod{ii} char(39) ...
                ' for the variable ' char(39) sMeta.metVar{ii} char(39) ...
                char(10) '.  The automated check believes these may be inconsistent.']);
        end
    end
    
end