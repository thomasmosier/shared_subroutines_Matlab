function shpOut = shp_MMap2Matlab(shpIn, varargin)

varMMap = 'ncst';
varDbf = 'dbf';
varDbfData = 'dbfdata';

varLon = 'X';
varLat = 'Y';


nShp = numel(shpIn.(varMMap)(:,1));

%Initialize struct
shpOut = struct;
%Preallocate:
shpOut(nShp).(varLon) = [];
shpOut(nShp).(varLat) = [];

for ii = 1 : nShp
    shpOut(ii).(varLon) = shpIn.(varMMap){ii}(:,2);
    shpOut(ii).(varLat) = shpIn.(varMMap){ii}(:,1);
end

%Variable input argument used to define fields to grab from input shapefile
if ~isempty(varargin(:))
    fldsShpIn = fieldnames(shpIn.(varDbf));
    for jj = 1 : numel(varargin(:))
        if iscell(varargin{jj})
            fldIn = varargin{jj}{1};
            fldOut = varargin{jj}{2};
        elseif ischar(varargin{jj})
            fldIn = varargin{jj};
            fldOut = fldIn;
        else
            warning('shpMMap2Matlab:unknownFldType', ['The variable input ' ...
                'argument should be char or cell. Class ' class(varargin{jj}) ...
                ' is not recognized.']);
            continue
        end
        indFld = find(strcmpi(fldsShpIn, fldIn));
        if ~isempty(indFld)
            shpOut(nShp).(fldOut) = [];
            
            for ii = 1 : nShp
                shpOut(ii).(fldOut) = shpIn.(varDbfData){ii,indFld};
            end
            clear ii
        else
            warning('shpMMap2Matlab:noField',['The field ' varargin{jj} ' was not found in the current shapefile.'])
        end
    end
    clear jj
end


