% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function strRes = t_res_geodata(time,timeUnits)

if numel(time) == 1
    strRes = 'unknown';
else
    dTime = mean(diff(time));

    if regexpbl(timeUnits,'day')
        if dTime > 20 && dTime < 35
            strRes = 'month';
        elseif dTime > 0.9 && dTime < 1.1
            strRes = 'day';
        elseif dTime > 359 && dTime < 366
            strRes = 'year';
        else
            strRes = 'unknown';
            warning('t_res_geodata:timeStep',['The average time step is ' num2str(dTime) ' does not correspond to a specific time step.']); 
        end
    elseif regexpbl(timeUnits,'hour')
        if dTime > 22 && dTime < 26
            strRes = 'day';
        elseif dTime < 22
            strRes = 'hour';
        elseif dTime > 600 && dTime < 750
            strRes = 'month';
        elseif dTime > 8616 && dTime < 8784
            strRes = 'year';
        else
            strRes = 'unknown';
            warning('t_res_geodata:timeStep',['The average time step is ' num2str(dTime) ' does not correspond to a specific time step.']); 
        end
    else
        strRes = 'unknown';
        warning('t_res_geodata:unknownUnits',['The units ' char(39) timeUnits char(39) ' are not recognized.']); 
    end
end

    
    