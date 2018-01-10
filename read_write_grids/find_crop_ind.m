function [indLonGp, indLatGp] = find_crop_ind(lonGrid, latGrid, lonBnds, latBnds, fr, type)


%Ensure lon and lat bounds are in correct order:
lonMin = min(lonBnds(:));
lonMax = max(lonBnds(:));
latMin = min(latBnds(:));
latMax = max(latBnds(:));

if strcmpi(type, 'in')
    indLonGp = [find(lonMin <= lonGrid, 1, 'first'), find(lonGrid <= lonMax, 1, 'last')];
    indLatGp = [find(latGrid <= latMax, 1, 'first'), find(latMin <= latGrid, 1, 'last')];
elseif strcmpi(type, 'out')
    indLonGp = [find(lonMin >= lonGrid, 1, 'last'), find(lonGrid >= lonMax, 1, 'first')];
    indLatGp = [find(latGrid >= latMax, 1, 'last'), find(latMin >= latGrid, 1, 'first')];
else
    error('cro_geodata:unknowntype', ['Method type ' type ' has not been programmed for.']);
end


%Check that two indices returned:
if numel(indLonGp) ~= 2 && any(isnan(indLonGp))
    error('cro_geodata:badLonInd', 'The lon indices returned are not correct.');
end
if numel(indLatGp) ~= 2 && any(isnan(indLatGp))
    error('cro_geodata:badLatInd', 'The lat indices returned are not correct.');
end

%add frame:
indLonGp = indLonGp + [-fr, +fr];
indLatGp = indLatGp + [-fr, +fr];

%Make sure indices are within bounds of available data:
if indLonGp(1) < 1
    indLonGp(1) = 1;
    warning('crop_geostrucct:outLonW','The Western lognitude cropping indice is out of bounds and therefore being set to 1.');
end
if indLatGp(1) < 1
    indLatGp(1) = 1;
    warning('crop_geostrucct:outLatN','The Northern latitude cropping indice is out of bounds and therefore being set to 1.');
end
if indLonGp(2) > numel(lonGrid)
    indLonGp(2) = numel(lonGrid);
    warning('crop_geostrucct:outLonE','The Eastern lognitude cropping indice is out of bounds and therefore being set to 1.');
end
if indLatGp(2) > numel(latGrid)
    indLatGp(2) = numel(latGrid);
    warning('crop_geostrucct:outLonE','The Southern latitude cropping indice is out of bounds and therefore being set to 1.');
end