function [ensoAnnBl, ensoAnnSig, ensoMnthSig] = geodata_ENSO(sData, varTas, lonEnso, latEnso, anomCrit)

%Calculate ENSO PHASE

%See the following reference for different ENSO region defintions:
%Trenberth, K. E. (1997). The definition of el nino. Bulletin of the American Meteorological Society, 78(12), 2771?2777.





%Definition is excursions from the decadal mean averaged over the Nino3.4
%region (5S-5N and 120-170W)
% lonBnds = [-170, -120];
% latBnds = [-5, 5];

nEnsoBs = 30*12; %Use 30 year baseline
nEnsoMnth = 5; %Use 5 month anomaly

varLon = 'longitude';
varLat = 'latitude';
varDate = 'date';

flag2Mat = 0;
if isstruct(sData)
    flag2Mat = 1;
    sData = {sData};
elseif ~iscell(sData)
    error('geodataEnso:unknwonDataType', ['The input data are in format ' class(sData) ', which is not programmed for.'])
end
    
ensoMnthSig = cell([numel(sData(:)), 1]);
ensoAnnBl = ensoMnthSig;
ensoAnnSig = ensoMnthSig;


%%Calculate input quantities to determine ENSO phase:
%For stock simulations
for kk = 1 : numel(sData(:))
    %Ensure data are monthly:
    if isfield(sData{kk}, 'timestep')
        if ~regexpbl(sData{kk}.timestep, 'month')
            error('geodataEnso:notMonthly','The input data do not appear to be monthly, which is a requirement.'); 
        end
    else
        if mode(diff(sData{kk}.(varDate)(:,2))) ~= 1
           error('geodataEnso:notMonthly','The input data do not appear to be monthly, which is a requirement.'); 
        end
    end
    
    %Ensure data cropped to ENSO region:
    [indLonGp, indLatGp] = find_crop_ind(sData{kk}.(varLon), sData{kk}.(varLat), lonEnso, latEnso, 0, 'out');
    sData{kk}.(varTas) = sData{kk}.(varTas)(:, indLatGp(1):indLatGp(end), indLonGp(1):indLonGp(end));
    
    yrsEnsoProj = sort(unique(sData{kk}.(varDate)(:,1)));
    nEnsoYrs = numel(yrsEnsoProj);
    nEnsoTs = numel(sData{kk}.(varDate)(:,1));

    %Initialize vector for years
    ensoAnnBl{kk} = zeros([nEnsoYrs, 1], 'single');
    ensoAnnSig{kk} = nan([nEnsoYrs, 1], 'single');

    %Calculate spatially-average temperature for ENSO region:
    ensoRegTas = nan([nEnsoTs, 1], 'single');
    for ii = 1 : nEnsoTs
        ensoRegTas(ii) = mean2d(squeeze(sData{kk}.(varTas)(ii,:,:))); 
    end
    clear ii

    %30-year moving average (based on monthly values)
    ensoRegRunAvg = runmean(ensoRegTas, ceil(nEnsoBs/2 - 1), [], 'edge');

    %5-month moving anomaly:
    ensoMnthSig{kk} = runmean(ensoRegTas, ceil(nEnsoMnth/2 - 1), [], 'edge') ...
        - ensoRegRunAvg;

    %Annual Enso Strength:
    for ii = 1 : nEnsoYrs
        indCurr = find(sData{kk}.date(:,1) == yrsEnsoProj(ii));

        %Calculate annual average strength:
        ensoAnnSig{kk}(ii) = nanmean(ensoMnthSig{kk}(indCurr));
        
        posEnsoTemp = ensoMnthSig{kk}(indCurr) >=  anomCrit;
        negEnsoTemp = ensoMnthSig{kk}(indCurr) <= -anomCrit;

        %Positive excursions:
        %Start of runs:
        I = run_length(posEnsoTemp);
        %Max length of positive run:
        nPosRun = max(I(posEnsoTemp == 1));
        %Negative excursions:
        %Strart of runs:
        I = run_length(negEnsoTemp);
        %Max length of positive run:
        nNegRun = max(I(negEnsoTemp == 1));
        
        if isempty(nPosRun)
            nPosRun = 0;
        end
        if isempty(nNegRun)
            nNegRun = 0;
        end
        
        if nPosRun > 5 && nNegRun < 5 %Positive
            ensoAnnBl{kk}(ii) =  1;
        elseif nNegRun > 5 && nPosRun < 5 %Negative
            ensoAnnBl{kk}(ii) = -1;
        end
    end
    clear ii
end

if flag2Mat == 1
    ensoAnnBl   =   ensoAnnBl{1};
    ensoAnnSig  = ensoAnnSig{1};
    ensoMnthSig = ensoMnthSig{1};
end
