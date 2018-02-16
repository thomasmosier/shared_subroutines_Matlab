function tasDelta = tas_range(tasMax, tasMin)

szTas = size(tasMax);
tasDelta = nan(szTas, 'single');

for ii = 1 : szTas(1)
    if ii ~= szTas(1)
        tasDelta(ii,:,:) = squeeze(tasMax(ii,:,:)) - squeeze(min(tasMin([ii,ii+1], :, :), [], 1));
    else
        tasDelta(ii,:,:) = squeeze(tasMax(ii,:,:)) - squeeze(tasMin(ii,:,:));
    end
end