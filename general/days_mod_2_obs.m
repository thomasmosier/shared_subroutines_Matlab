function [nDaysMod, nDaysObs, indModStrt, indModEnd] = days_mod_2_obs(daysModStrt, daysModEnd, daysObsStrt, daysObsEnd)

nObs = numel(daysObsStrt(:,1));
nDaysMod   = nan(nObs, 1);
nDaysObs   = nDaysMod;
indModStrt = nDaysMod;
indModEnd  = nDaysMod;

%Loop over observations
for ii = 1 : nObs
    %Find days seperating mod and obs
    nDatesMod = numel(daysModStrt(:,1));
    daysStrtTemp = nan(nDatesMod,1);
    daysEndTemp = daysStrtTemp;
    for zz = 1 : nDatesMod
        daysStrtTemp(zz) = daysModStrt(zz) - daysObsStrt(ii);
        daysEndTemp(zz)  = daysObsEnd(ii) - daysModEnd(zz);
    end

    %Find indices of model that best include
    %observations (going over if necessary)
    %For start:
    if any(daysStrtTemp <= 0)
        indModStrt(ii) = find(daysStrtTemp <= 0, 1, 'last');
    else
        indModStrt(ii) = 1;
    end
    %For end:
    if any(daysEndTemp <= 0)
        indModEnd(ii) = find(daysEndTemp <= 0, 1, 'first');
    else
        indModEnd(ii) = numel(daysEndTemp);
    end

    %Interpolate model to observation dates:
    nDaysMod(ii) = daysModEnd(indModEnd(ii)) - daysModStrt(indModStrt(ii));
    nDaysObs(ii) = daysObsEnd(ii) - daysObsStrt(ii);
end