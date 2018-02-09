function blOut = align_check(lon1, lon2, lat1, lat2)

if ~any(size(lon1) == 1)
    error('alignChck:dimLon1', 'Neither of the lon1 dimensions is 1. This function is designed to work with vector inputs.');
elseif ~any(size(lon2) == 1)
    error('alignChck:dimLon2', 'Neither of the lon2 dimensions is 1. This function is designed to work with vector inputs.');
elseif ~any(size(lat1) == 1)
    error('alignChck:dimLat1', 'Neither of the lat1 dimensions is 1. This function is designed to work with vector inputs.');
elseif ~any(size(lat2) == 1)
    error('alignChck:dimLat2', 'Neither of the lat2 dimensions is 1. This function is designed to work with vector inputs.');
end
    
dLat = nanmean([abs(diff(lat1(:))); abs(diff(lat2(:)))]);
dLon = nanmean([abs(diff(lon1(:))); abs(diff(lon2(:)))]);

ordLat = -order(dLat)+1;
ordLon = -order(dLon)+1;


if ~isequal(round2(lon1, ordLon), round2(lon2, ordLon))
    blOut = 0;
elseif ~isequal(round2(lat1, ordLat), round2(lat2, ordLat))
    blOut = 0;
else
    blOut = 1;
end