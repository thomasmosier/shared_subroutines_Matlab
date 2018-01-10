function q0 = humidity_sat_specific(tas, ps, unitsTas, unitsPs)
%Calculate saturated specific humidity based on:

%tas = 2m air temperature (units = Celsius or Kelvin)
%ps = Surface air pressure (units = Pascals)
%q0 = unitless (kg / kg)

Tt = 273.16; %Triple point temperature (Kelvin)
eps = 0.62198;
Tk = 273.15;

if regexpbl(unitsTas, 'Celsius')
    tas = tas + Tk;
elseif ~regexpbl(unitsTas, 'kelvin') && ~strcmpi(unitsTas, 'K')
   error('humidity_sat_specific:UnknownTasUnits', ['The units of ' unitsTas ...
       ' have not been programmed for.']); 
end

if ~regexpbl(unitsPs, 'Pa')
    error('humidity_sat_specific:UnknownPsUnits', ['The units of ' unitsPs ...
       ' have not been programmed for.']); 
end

indRng = (1:numel(tas));
% indRng = find(223 < tas & tas < 373);

% indLow = find(173 < tas & tas < Tk);
% indHgh = find(Tk <= tas & tas < 373);

if nanmin(min2d(tas)) < 223 || nanmax(max2d(tas)) > 373
    warning('humidity_sat_specific:tasOutRang', ['Some temperature values '...
        'are outside the acceptable range (223 to 373). Min = ' ...
        num2str(nanmin(min2d(tas))) '; Max = ' num2str(nanmax(max2d(tas)))]); 
end

logEs = nan(size(tas), 'single');

logEs(indRng) = 10.79574*(1-Tt./tas) - 5.028*log10(tas/Tt) ...
    + 1.50475*10^(-4)*(1-10.^(-8.2969*(tas/Tt - 1))) ...
    + 0.42873*10^(-3)*(10.^(4.76955*(1-Tt./tas)) - 1) ...
    + 0.78614 + 2.0;

%CODE TO TEST DIFFERENCE BASED ON FORMULA:
% %Formulation based on Clasius-Clapyron equation:
% a = 17.67; %unitless
% b = 243.5; %Celsius
% esCC = 100*6.122*exp(a*(tas - Tk)./(b+(tas - Tk))); %units = Pascals
% diff = (esCC - 10.^(logEs)) ./ ((esCC + 10.^(logEs))/2);
%   
% disp(['Mean difference: ' num2str(mean(mean2d(diff))) ...
%     '; Max difference: ' num2str(max(max2d(diff)))]);


%saturation vapor pressure over ice:
% logEs(indHgh) = -9.09685*(Tt/tas - 1) - 3.56654*log10(Tt/tas) ...
%     + 0.87682*(1 - tas/Tt) + 0.78614 + 2.0;

%Calculate specific saturation humidity
q0 = eps * 10.^(logEs) ./ ps;