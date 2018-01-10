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

function varargout = same_grid(lonPr,latPr, type, varargin)

varargout = varargin(:);

for ii = 1 : numel(varargin(:))
    if ~isequal(round2(lonPr,4), round2(varargin{ii}.lon,4)) || ~isequal(round2(latPr,4), round2(varargin{ii}.lat,4))
        if numel(size(varargin{ii}.data)) == 2
            dataRe = nan(numel(latPr),numel(lonPr),'single');
        elseif numel(size(varargin{ii}.data)) == 3
            dataRe = nan(numel(varargin{ii}.time),numel(latPr),numel(lonPr),'single');
        else
            
        end

        
        if regexpbl(type,'area')
            sRef.lat = latPr;
            sRef.lon = lonPr;
            sData = regrid_geodata(varargin{ii},sRef);
            varargout{ii}.data = sData.dataRe;
            varargout{ii}.lon = sData.lonRe;
            varargout{ii}.lat = sData.latRe;
            clear sData
        else
            for mm = 1 : numel(varargin{ii}.time)
                if regexpbl(type,'pchip')
                    if numel(size(varargin{ii}.data)) == 2
                        dataRe = PCHIP_2D(varargin{ii}.lon, varargin{ii}.lat, varargin{ii}.data, lonPr,latPr);
                    elseif numel(size(varargin{ii}.data)) == 3
                        dataRe(mm,:,:) = PCHIP_2D(varargin{ii}.lon, varargin{ii}.lat, squeeze(varargin{ii}.data(mm,:,:)), lonPr,latPr);
                    else
                        error('same_grid:dim','Arrays with 2 or 3 dimensions are expected');
                    end
                else
                    if numel(size(varargin{ii}.data)) == 2
                        dataRe = interp2(varargin{ii}.lon, varargin{ii}.lat, varargin{ii}.data, lonPr,latPr,type);
                    elseif numel(size(varargin{ii}.data)) == 3
                        dataRe(mm,:,:) = interp2(varargin{ii}.lon, varargin{ii}.lat, squeeze(varargin{ii}.data(mm,:,:)), lonPr,latPr,type);
                    else
                        error('same_grid:dim','Arrays with 2 or 3 dimensions are expected');
                    end
                    
                end
            end
            
            varargout{ii}.data = dataRe;
            varargout{ii}.lon = lonPr;
            varargout{ii}.lat = latPr;
        end
    end
end
