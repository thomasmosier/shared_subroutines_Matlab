function [blOut, lonO, latO, varargout] = crd_within(lon1, lon2, lat1, lat2, varargin)

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


nLonSame = numel(intersect(round2(lon1(:), ordLon), round2(lon2(:), ordLon)));
nLatSame = numel(intersect(round2(lat1(:), ordLat), round2(lat2(:), ordLat)));

if nLatSame >= min([numel(lat1), numel(lat2)]) && nLonSame >= min([numel(lon1), numel(lon2)]) 
	blOut = 1; %One grid contained in the other
else
	blOut = 0; %Grid mostly disimilar
end


latComb = [lat1(:) ; lat2(:) ];
lonComb = [lon1(:)', lon2(:)'];

[~, iLatO, ~] = unique(round2(latComb, ordLat));
[~, iLonO, ~] = unique(round2(lonComb, ordLon));

latO = latComb(iLatO);
lonO = lonComb(iLonO);
    

if blOut == 1 && ~isempty(varargin(:))
    data1 = varargin{1};
    data2 = varargin{2};

        
    if numel(lat1) ~= numel(lat2)
        if numel(lat1) < numel(lat2)
            tempO = nan([numel(data1(:,1,1)), numel(lat2), numel(data1(1,1,:))]);
            [~, indLat2] = ismember(round2(lat1, ordLat), round2(lat2, ordLat));
            tempO(:,indLat2,:) = data1;
            data1 = tempO;
        else
            tempO = nan([numel(data2(:,1,1)), numel(lat1), numel(data2(1,1,:))]);
            [~, indLat1] = ismember(round2(lat2, ordLat), round2(lat1, ordLat));
            tempO(:,indLat1,:) = data2;
            data2 = tempO;
        end
    end
    
    
    if numel(lon1) ~= numel(lon2)
        if numel(lon1) < numel(lon2)
            tempO = nan([numel(data1(:,1,1)), numel(data1(1,:,1)), numel(lon2)]);
            [~, indLon2] = ismember(round2(lon1, ordLon), round2(lon2, ordLon));
            tempO(:,:,indLon2) = data1;
            data1 = tempO;
        else
            tempO = nan([numel(data2(:,1,1)), numel(data2(1,:,1)), numel(lon1)]);
            [~, indLon1] = ismember(round2(lon2, ordLon), round2(lon1, ordLon));
            tempO(:,:,indLon1) = data2;
            data2 = tempO;
        end
    end
    
    if nargout == 5
        varargout{1} = data1;
        varargout{2} = data2;
    end
end