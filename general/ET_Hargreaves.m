function ET = ET_Hargreaves(tMax, tMin, lat, date, Kt)

% %Solar constant
% I0 = 1360; %(W/m^2) DeWalle and Rango, Principles of Snow Hydrology
% 
% raad = nan(numel(date(:,1)), 1);
% for ii = 1 : numel(date(:,1))
%     %Julian day
%     dayJ = Julian_day(date,'ofyear');
% 
%     %Potential sunlight hours (assume flat)
%     angleDec = 23.45*sind((dayJ+284)*360/365);
% 
%     hrs = real(2*acosd(-tand(lat)*tand(angleDec(ii)))/15);
%         hrs(hrs > 24) = 24;
%         hrs(hrs <  0) =  0;
%         hrs = squeeze(hrs);
% 
%     %Eccentricity of Earth's orbit:
%     %ecc = (1 - e^2)/(1+e*cosd(theta))
%     ecc = 1 + 0.033*cosd(dayJ*360/365);
% 
%     hourAngle1 = -7.5*hrs;
%     hourAngle2 = 7.5*hrs;
% 
%     %Integrated cos(Z)
%     cZenith = hrs.*sind(lat).*sind(angleDec(ii)) ...
%         + (1/15)*cosd(lat).*cosd(angleDec(ii)).*(sind(hourAngle2)-sind(hourAngle1)); %Unitless
%     %Daily Radiation (Units of Watts / meter squared)
%     raad(ii) = (I0./ecc)*cZenith/24;
%     
% end
% 
% %Convert raadiation to MJ/day
% facW2MJ = (24*60*60)/(10^6);
% raadMJDy = facW2MJ*raad;
% 
% facMJ2MM = 0.408;
% raadMMDy = facMJ2MM*raadMJDy;

raadLat = [38, 36, 34];
raadTable = [6.9, 9.0, 11.8, 14.5, 16.4, 17.2, 16.7, 15.3, 12.8, 10.0, 7.5, 6.1; ...
             7.4, 9.4, 12.1, 14.7, 16.4, 17.2, 16.7, 15.4, 13.1, 10.6, 8.0, 6.6; ...
             7.9, 9.8, 12.4, 14.8, 16.5, 17.1, 16.8, 15.5, 13.4, 10.8, 8.5, 7.2];

if lat > max(raadLat) || lat < min(raadLat)
    error('etHargreaves:latRange',['Latitude ' num2str(lat) ' is outside the range of programmed latitudes (' num2str(min(raadLat)) ' to ' num2str(max(raadLat)) ')' ]);
end

[~, indLat] = min(abs(raadLat - lat));
raadMMDy = nan(numel(date(:,1)), 1);         
for ii = 1 : numel(date(:,1))
    raadMMDy(ii) = raadTable(indLat, date(ii,2));
end

tRng = [tMax(:), tMin(:)];
tMaxUse = max(tRng, [], 2);
tMinUse = min(tRng, [], 2);

tAvg = mean([tMaxUse, tMinUse], 2);

%Hargreave's formula:
ET = 0.0135*Kt*(tAvg+17.78).*raadMMDy.*(tMaxUse-tMinUse).^(0.5);