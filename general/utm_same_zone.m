function [xOut, yOut] = utm_same_zone(xIn, yIn, zoneIn, ZoneForce)

%Check sizes of inputs match
if ~all(size(xIn) == size(yIn)) || ~all(size(xIn) == size(zoneIn))
    error('utm_same_zone:inputSizeDiff', 'The sizes of the input arrays do not match.');
end

%Initialize output:
xOut = nan(size(xIn));
yOut = nan(size(yIn));

%Loop over input points:
for ii = 1 : numel(xIn)
    [lat, lon] = utm2deg(xIn(ii), yIn(ii), zoneIn{ii});
    
    [xOut(ii), yOut(ii)] = deg2utm_force_zone(lat, lon, char(ZoneForce));
end