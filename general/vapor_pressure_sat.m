function eSat = vapor_pressure_sat(tas, unitsTas)
%Calculate saturated specific humidity based on:

%tas = 2m air temperature (units = Celsius or Kelvin)
%ps = Surface air pressure (units = Pascals)
%q0 = unitless (kg / kg)

Tk = 273.15;


if regexpbl(unitsTas, 'kelvin') || strcmpi(unitsTas, 'K')
    tas = tas - Tk;
elseif ~regexpbl(unitsTas, 'Celsius')
   error('humidity_sat_specific:UnknownTasUnits', ['The units of ' unitsTas ...
       ' have not been programmed for.']); 
end

a = 17.67; %unitless
b = 243.5; %Celsius

eSat = 100*6.122*exp(a*tas./(b+tas)); %units = Pascals
