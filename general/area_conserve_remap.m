function dataO = area_conserve_remap(lonI, latI, dataI, lonO, latO, varargin)


%Default, set bins = 1
bin = 1;
if ~isempty(varargin(:))
   for ii = 1 : numel(varargin(:))
      if regexpbl(varargin{ii}, 'bin')
          bin = varargin{ii+1};
      end
   end
end

if ~isequal(bin, round2(bin, 2))
   warning('areaConserveRemap:binVal', ['Only integer values of the '...
       'bin parameter are accepted. Therefore ' num2str(bin) ' is being '...
       'rounded to ' num2str(round(bin)) '.' ]);
   bin = round(bin);
end

%Area parameters (used in aggregation):
rEarth = 6371000;
const = (pi/180)*rEarth^2;

%Compute box "interfaces" / edges:
fLonI = box_edg(lonI);
fLatI = box_edg(latI);

fLonO = box_edg(lonO);
fLatO = box_edg(latO);

edgLonIW = fLonI(1:end-1);
edgLonIE = fLonI(2:end);

edgLatIN = fLatI(1:end-1);
edgLatIS = fLatI(2:end);


%Calculate edges of output grid using bin factor:
if bin == 1
    edgLonOW = fLonO(1:end-1);
    edgLonOE = fLonO(2:end);

    edgLatON = fLatO(1:end-1);
    edgLatOS = fLatO(2:end);
else
    edgLonOW = nan(bin*numel(fLonO(2:end)),1);
    edgLonOE = edgLonOW;
    
    cntr = 1;
    for ii = 1 : numel(fLonO(2:end))
        temp = linspace(fLonO(ii),fLonO(ii+1),bin+1);
        edgLonOW(cntr:cntr+bin-1) = temp(1:end-1);
        edgLonOE(cntr:cntr+bin-1) = temp(2:end);
        
        cntr = cntr + bin;
    end
    clear ii
    
    
    edgLatON = nan(bin*numel(fLatO(2:end)),1);
    edgLatOS = edgLatON;
    
    cntr = 1;
    for ii = 1 : numel(fLatO(2:end))
        temp = linspace(fLatO(ii),fLatO(ii+1),bin+1);
        edgLatON(cntr:cntr+bin-1) = temp(1:end-1);
        edgLatOS(cntr:cntr+bin-1) = temp(2:end);
        
        cntr = cntr + bin;
    end
    clear ii    
end


%Make LON vector mapping from input to output:
ixI = zeros(4*max([numel(lonI),numel(lonO)]),1, 'uint16');
ixO = ixI;
dx  = nan(size(ixI), 'single');
xW = nan(size(ixI), 'single');
xE = xW;
cntrX = 0;
for i = 1 : numel(lonO)
    for ii = 1 : numel(lonI)
        if edgLonIE(ii) > edgLonOW(i) && edgLonOE(i) > edgLonIW(ii)
            cntrX = cntrX + 1;
            
            ixI(cntrX) = ii;
            ixO(cntrX) =  i;
            [dx( cntrX), indX] =  min([ ...
                edgLonIE(ii)-edgLonIW(ii), ...
                edgLonOE( i)-edgLonOW( i), ...
                edgLonIE(ii)-edgLonOW( i), ...
                edgLonOE( i)-edgLonIW(ii)]);
            
            switch indX
                case 1
                    xW(cntrX) = edgLonIE(ii);
                    xE(cntrX) = edgLonIW(ii);
                case 2
                    xW(cntrX) = edgLonOE( i);
                    xE(cntrX) = edgLonOW( i);
                case 3
                    xW(cntrX) = edgLonIE(ii);
                    xE(cntrX) = edgLonOW( i);
                case 4
                    xW(cntrX) = edgLonOE( i);
                    xE(cntrX) = edgLonIW(ii);
            end
        end
    end
end
ixI = ixI(1:cntrX);
ixO = ixO(1:cntrX);
dx  =  dx(1:cntrX);
xW  =  xW(1:cntrX);
xE  =  xE(1:cntrX);

%Make LAT vector mapping from input to output:
iyI = zeros(4*max([numel(latI),numel(latO)]),1, 'uint16');
iyO = iyI;
dy  = nan(size(iyI), 'single');
yN = nan(size(iyI), 'single');
yS = yN;
cntrY = 0;
for j = 1 : numel(latO)
    for jj = 1 : numel(latI)
        if edgLatIN(jj) > edgLatOS(j) && edgLatON(j) > edgLatIS(jj)
            cntrY = cntrY + 1;
            
            iyI(cntrY) = jj;
            iyO(cntrY) =  j;
            [dy( cntrY), indY] =  min(abs([...
                edgLatON(j)-edgLatOS(j), ...
                edgLatIN(jj)-edgLatIS(jj), ...
                edgLatIN(jj)-edgLatOS(j), ...
                edgLatON(j)-edgLatIS(jj)]));
            
            switch indY 
                case 1
                    yN(cntrY) = edgLatON(j);
                    yS(cntrY) = edgLatOS(j);
                case 2
                    yN(cntrY) = edgLatIN(jj);
                    yS(cntrY) = edgLatIS(jj);
                case 3
                    yN(cntrY) = edgLatIN(jj);
                    yS(cntrY) = edgLatOS(j);
                case 4
                    yN(cntrY) = edgLatON(j);
                    yS(cntrY) = edgLatIS(jj);
            end
        end
    end
    clear jj
end
clear j
iyI = iyI(1:cntrY);
iyO = iyO(1:cntrY);
dy  =  dy(1:cntrY);
yN  =  yN(1:cntrY);
yS  =  yS(1:cntrY);


%Aggregate:
%area formula = (pi/180)*rEarth^2*abs(sind(latBN)-sind(latBS)).*abs(lonBE-lonBW);
dataO = zeros(numel(latO), numel(lonO), 'single');
for j = 1 : cntrY
    indYO = iyO(j);
    indYI = iyI(j);
    areaY = abs(sind(yN(j))-sind(yS(j)));
    
    for i = 1 : cntrX
        indXO = ixO(i);
        indXI = ixI(i);
        
        dataO(indYO, indXO) = ...
              dataO(indYO, indXO) ...
            + dataI(indYI, indXI)*(const*areaY.*abs(xW(i)-xE(i)));   
%         dataO(indYO, indXO) = ...
%               dataO(indYO, indXO) ...
%             + dataI(indYI, indXI)*dy(j)*dx(i);
    end
    clear i
end
clear j

%Normalize:
areaO = area_geodata(fLonO,fLatO, 'e');
dataO = dataO ./ areaO;

%Set missing values to nan:
[meshX, meshY] = meshgrid(ixO, iyO);
indNan = setdiff((1:numel(dataO)), unique(sub2ind(size(dataO), meshY(:), meshX(:))));

if ~isempty(indNan)
    dataO(indNan) = nan;
end
