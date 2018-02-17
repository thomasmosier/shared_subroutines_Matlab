function tasRange = tas_range(tasMax, tasMin)

szTas = size(tasMax);
tasRange = nan(szTas, 'single');

for ii = 1 : szTas(1)
    if ii ~= szTas(1)
        tasRange(ii,:,:) = squeeze(tasMax(ii,:,:)) - squeeze(min(tasMin([ii,ii+1], :, :), [], 1));
    else
        tasRange(ii,:,:) = squeeze(tasMax(ii,:,:)) - squeeze(tasMin(ii,:,:));
    end
end