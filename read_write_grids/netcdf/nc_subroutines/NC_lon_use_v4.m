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

function [lonTr, indLonUse] = NC_lon_use_v4(lonOrg, attLon, sMeta)

%Make sure col vector
lonOrg = lonOrg(:)';

%%Find indices for longitude of interest in NetCDF or GRB file:
indLonUnt = strcmpi(attLon,'units');
[rowLon, colLon] = find(indLonUnt == 1);

%Checking units may or may not be important:
lonUnits = attLon{rowLon, colLon+1};

%Determine frame to use (around coordinates bounds)
if isfield(sMeta,'frame') 
    frame = sMeta.frame;
else
    frame = 0; 
end

if isfield(sMeta,'crop') 
    cropType = sMeta.crop;
elseif isfield(sMeta,'cropType') 
    cropType = sMeta.cropType;
else
    cropType = 'out'; 
end

%Find longitude indices that are within the desired region:
%case where 'degrees East' are used:
if ~isempty(regexpi(lonUnits,'deg')) && ~isempty(regexpi(lonUnits,'east'))
    if max(lonOrg) > 180 && min(lonOrg) >= 0 %Longitude goes from 0 to 360 degrees east; otherwise assume -180 -> 180
        %Transform to -180 -> 0 -> 180
        indOrg180 = find(180 <= lonOrg,1, 'first');
        lonTr = [lonOrg(indOrg180:end)-360, lonOrg(1:indOrg180-1)];

        if isfield(sMeta,'crd') && all(~isnan(sMeta.crd(1:2)))
            lonMin = min(sMeta.crd(1:2));
            lonMax = max(sMeta.crd(1:2));
            
            lonTrEdg = box_edg(lonTr);
            
            if strcmpi(cropType, 'in')
                indW = find(lonTrEdg(1:end-1) >= lonMin , 1, 'first');
                indE = find(lonMax >= lonTrEdg(2:end), 1, 'last');
            elseif strcmpi(cropType, 'out')
                indW = find(lonMin >= lonTrEdg(1:end-1), 1, 'last');
                indE = find(lonTrEdg(2:end) >= lonMax, 1, 'first');
            else
                error('cro_geodata:unknowntype', ['Method type ' cropType ' has not been programmed for.']);
            end
            
%             lonTrEdg(indW)
%             lonTrEdg(indE+1)
            
            if isempty(indW)
                indW = 1; 
            end
            if isempty(indE)
                indE = numel(lonOrg); 
            end
            
            indW = indW - frame;
            indE = indE + frame;
            
            if indW < 1
                indW = 1;
                warning('crop_geostrucct:outLonW','The Northern latitude cropping indice is out of bounds and therefore being set to 1.');
            end
            if indE > numel(lonOrg)
                indE = numel(lonOrg);
                warning('crop_geostrucct:outLonE','The Northern latitude cropping indice is out of bounds and therefore being set to 1.');
            end
            
%             indLonUse = (indW : indE);
            
%             %Find indices corresponding to bounds of selection (in transformed 
%             %coordinates):
%             indW = find(lonTrEdg(1:end-1) < sMeta.crd(1), 1, 'last');
%             if isempty(indW) && 2*lonTrEdg(1) - lonTrEdg(2) <= -179 && 2*lonTrEdg(end) - lonTrEdg(end-1) > 179   %Data extends beyond -180
%                indW = find(lonTrEdg(1:end-1) -360 < sMeta.crd(1), 1, 'last') - numel(lonOrg) - 1;
%             elseif isempty(indW)
%                 indW = 1;
%                 flagLon(1) = 1;
%             end
%             indE = find(lonTrEdg(2:end) > sMeta.crd(2), 1, 'first');
%             if isempty(indE) && 2*lonTrEdg(1) - lonTrEdg(2) <= -179 && 2*lonTrEdg(end) - lonTrEdg(end-1) > 179 %Data extends beyond 180
%                indE = find(lonTrEdg(2:end) + 360 > sMeta.crd(2), 1, 'first') + numel(lonOrg);
%             elseif isempty(indE)
%                 indE = numel(lonOrg);
%                 flagLon(2) = 1;
%             end

            %Identify indice where transformed coordinates are closest to 0:
            indTr0 = find(lonTr >= 0, 1, 'first');
%             lonTr(indTr0)
            
            %Rearange indices based on GCM (original) orientation:
            if indW < indTr0 && indE < indTr0 %Both negative (Western hemisphere)
                indLonUse = (numel(lonOrg) - indOrg180 + indW + 1: numel(lonOrg) - indOrg180 + indE + 1); %Empirically checked
            elseif indW < indTr0 && indE >= indTr0 %western neg and eastern positive
                indLonUse = [(numel(lonOrg) - indOrg180 + indW + 1 : numel(lonOrg)), (1 : indE - indOrg180 + 1)]; %Empirically checked
            elseif indW >= indTr0 && indE >= indTr0 %Both positive (Eastern hemisphere)
                indLonUse = (indW - indOrg180 + 1 : indE - indOrg180 + 1); %Empirically checked
            else
                error('NC_lon_use:lonOrientation','Case unexpected and not coded for.');
            end

            %If longitude on edge of grid being used, must wrap:
            indLonUse = bnd_wrap(indLonUse, numel(lonOrg));

            %Select longitude values for output:
%             lonTr = lonTr(indLonUse);
            lonTr = lonOrg(indLonUse);
        else
            indLonUse = (1:numel(lonOrg));
        end
        
        lonTr( lonTr >= 180 ) = lonTr( lonTr >= 180 ) - 360;
    else %Longitude goes from -180 to 180 degrees east
        if isfield(sMeta,'crd')  && all(~isnan(sMeta.crd(1:2)))
            lonMin = min(sMeta.crd(1:2));
            lonMax = max(sMeta.crd(1:2));
            
            lonOrgEdg  = box_edg(lonOrg);

            if strcmpi(cropType, 'in')
                indW = find(lonMin <= lonOrgEdg(1:end-1), 1, 'first');
                indE = find(lonOrgEdg(2:end) <= lonMax, 1, 'last');
            elseif strcmpi(cropType, 'out')
                indW = find(lonMin >= lonOrgEdg(1:end-1), 1, 'last');
                indE = find(lonOrgEdg(2:end) >= lonMax, 1, 'first');
            else
                error('cro_geodata:unknowntype', ['Method type ' cropType ' has not been programmed for.']);
            end
            
            if isempty(indW)
                indW = 1; 
            end
            if isempty(indE)
                indE = numel(lonOrg); 
            end
            
            indW = indW - frame;
            indE = indE + frame;
            
            if indW < 1
                if indW + frame < 1
                   warning('crop_geostrucct:outLonW', ['The Western ' ...
                       'longitude cropping indice is out of bounds and is therefore being set to 1.']); 
                end
                indW = 1;
            end
            if indE > numel(lonOrg)
                if indE - frame > numel(lonOrg)
                    warning('crop_geostrucct:outLonE',['The Eastern '...
                        'longitude cropping indice is out of bounds and '...
                        'is therefore being set to ' num2str(numel(lonOrg)) '.']);
                end
                indE = numel(lonOrg);
            end
            
            indLonUse = (indW : indE);
            
%             indW = find(lonOrgEdg(1:end-1) <= sMeta.crd(1), 1, 'last');
%             if isempty(indW) && abs(diff([lonOrgEdg(1),sMeta.crd(1)])) < 0.1*nanmean(abs(diff(lonOrgEdg)))
%                 indW = 1;
%             elseif isempty(indW) && 2*lonOrgEdg(1) - lonOrgEdg(2) <= -179 && 2*lonOrgEdg(end) - lonOrgEdg(end-1) > 179   %Data extends beyond -180
%                 indW = find(lonOrgEdg(1:end-1) -360 < sMeta.crd(1), 1, 'last') - numel(lonOrg) - 1;
%             elseif isempty(indW)
%                 indW = 1;
%                 flagLon(1) = 1;
%             end
%             indE = find(lonOrgEdg(2:end) > sMeta.crd(2), 1, 'first');
%             if isempty(indE) && abs(diff([lonOrgEdg(end),sMeta.crd(2)])) < 0.1*nanmean(abs(diff(lonOrgEdg)))
%                 indE = 1;
%             elseif isempty(indE) && 2*lonOrgEdg(1) - lonOrgEdg(2) <= -179 && 2*lonOrgEdg(end) - lonOrgEdg(end-1) > 179 %Data extends beyond 180
%                indE = find(lonOrgEdg(2:end) + 360 > sMeta.crd(2), 1, 'first') + numel(lonOrg); 
%             elseif isempty(indE)
%                 indE = numel(lonOrg);
%                 flagLon(2) = 1;
%             end

%             if ~any(flagLon)
%                 indLonUse = (indW - frame : indE + frame);
%             elseif flagLon(1) && ~flagLon(2)
%                 indLonUse = (indW : indE + frame);
%             elseif ~flagLon(1) && flagLon(2)
%                  indLonUse = (indW - frame : indE);
%             else
%                 indLonUse = (indW : indE);
%             end

            %If longitude on edge of grid being used, must wrap:
            indLonUse = bnd_wrap(indLonUse, numel(lonOrg));
        else
            indLonUse = (1:numel(lonOrg));
        end
        %Select longitude values for output:
        lonTr = lonOrg(indLonUse);
%     else
%         keyboard
%         error('gcm_load:coordOffset',['An unknown longitudinal offset '...
%             'is used in the current GCM time series.' char(10) ...
%             'There is likely nothing wrong except that this case must be coded for.']);
    end
    
    %Bnd Wrap may introduce discontinuity, so remove unnatural jumps
    indRem = find(abs(diff(lonTr))/(5*mean(abs(diff(lonTr)))) > 1);
    if ~isempty(indRem)
        warning('NC_lon_use:discontinuity','Output grid may have longitudinal discontinuity. If so, please remove.');
%         if diff([lonOrg(indLonUse(indRem(1))),lonOrg(indLonUse(indRem(1)+1))]) < 0
%             if indRem(1) > 0.5*numel(indLonUse)
%                 indLonUse = indLonUse(1:indRem(1));
%             else
%                 indLonUse = indLonUse(indRem(1):end);
%             end
%         end
    end

    %If lonTr wraps backwards across -180, adjust values so they extend
    %beyond -180 (necessary for interpolation on flat surface)
    indEeg = find(diff(lonTr) < 0);
    if ~isempty(indEeg)
        if numel(indEeg) == 1
            if indEeg < 0.5*numel(lonTr)
                lonTr(1:indEeg) = lonTr(1:indEeg) - 360;
            else
                lonTr(indEeg+1:end) = lonTr(indEeg+1:end) + 360;
            end      
        elseif numel(indEeg) == 2
            lonTr(1:indEeg(1)) = lonTr(1:indEeg(1)) - 360;
            lonTr(indEeg(2)+1:end) = lonTr(indEeg(2)+1:end) + 360;
        else
            error('NC_lon_use:wrap2Ind','This case not coded for.')
        end
    end
else
    error('gcm_load:coordSys','An unknown longitudinal coordinate system is used in the current GCM time series.');
end