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

function [varargout] = days_2_date_v2(daysIn, dateRef, strCal)



%Date Vectors have the form: [year, month, day].  
%Can handle multiple dates (where each date is it's own row)

%Change strCal into numerical value defined in 'cal_class'
blCal = cal_class(strCal);

if any(isnan(dateRef))
    warning('days2Date:unknwonRef','Function returning early because reference date contains nan values.');
    varargout = nan([numel(daysIn), numel(dateRef)]);
    return
end
    
%Make row vector if column:
if length(dateRef(:,1)) > 1 && length(dateRef(1,:)) == 1
    dateRef = dateRef';
end
if length(daysIn(1,:)) > 1 && length(daysIn(:,1)) == 1
    daysIn = daysIn';
end

%Check that 'daysIn' is vector and not matrix
if sum(size(daysIn)>1) >1
   error('days_2_date:daysMatrix',[char(39) 'DaysIn' char(39) ' must be vector, not matrix.']);
end

% if any(daysIn < 0)
%    error('days2Date:negInput','Input days must be positive.') 
% end
    
%Initialize date output vector:
nPrec = numel(dateRef(1,:));
dateOut = nan(numel(daysIn), nPrec, 'single');

if blCal == 0 || blCal == 1 || blCal == 3 %Gregorian, 365 (no leap), or Julian
    daysMnths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    daysYr = 365;
elseif blCal == 2 %360 days
    daysMnths = 30*ones(1,12);
    daysYr = 360;
else
    error('days2Date:unknownCal',['The calendar type ' strCal ' is not recognized.']);
end
daysMnthCum = cumsum(daysMnths);
daysMnthCumLp = [daysMnthCum(1), daysMnthCum(2:end) + 1];

%Modify "days" so that reference is 1st of year
switch blCal 
    case 0 %Gregorian
        dateRefTemp = ones(1,6);
        dateRefTemp(1:numel(dateRef)) = dateRef;

        dateRefJan1 = ones(size(dateRefTemp));
        dateRefJan1(1) = dateRef(1);
        daysExtra = datenum(datestr(dateRefTemp)) - datenum(datestr(dateRefJan1));
    case 1 %No leap
        daysExtra = sum(daysMnths(1:dateRef(2)-1)) + dateRef(3)-1;
    case 3 %Julian
        if mod(dateRef(1), 4) && dateRef(2) > 2
            daysExtra = sum(daysMnths(1:dateRef(2)-1)) + dateRef(3);  
        else
            daysExtra = sum(daysMnths(1:dateRef(2)-1)) + dateRef(3)-1;
        end
    case 2 %360
        daysExtra = (dateRef(2)-1)*30 + dateRef(3)-1;
end


daysIn = daysIn + daysExtra;

if ~any(isnan(dateRef))
    for ii = 1 : numel(daysIn)
        %If 'daysIn' is NaN, output = nan
        if isnan(daysIn(ii))
            continue
        end

        
        %This is a workaround for error in leap day encoding:
        if ii == 28 && numel(daysIn) == 29 && all(diff(daysIn(1:27)) == 1) && diff(daysIn(28:29)) > 10000
            %nYrs = nYrs (same as last value)
            daysIn(29) = daysIn(28) + 1;
            nYrs = floor(daysIn(ii)/daysYr);
        else
            %Find year no counting leap:
            nYrs = floor(daysIn(ii)/daysYr);
        end
        
        %Find leaps
        switch blCal 
            case 0 %Gregorian
                lps = is_leap_v2(dateRef(1) + (0:nYrs));
            case 3 %Julian
                lps = mod(dateRef(1) + (0:nYrs), 4);
            otherwise
                lps = zeros(2,1);
        end

        %See if number of leaps changes number of years:
        if nYrs > 0
            nLp = sum(lps(1:end-1));
            
            if nYrs*daysYr+nLp > daysIn(ii)
                nYrs = nYrs - 1;

                lpCurr = lps(end-1);

                daysIn(ii) = daysIn(ii) - nYrs*daysYr - nLp + lpCurr;
            else
                lpCurr = lps(end);

                daysIn(ii) = daysIn(ii) - nYrs*daysYr - nLp;
            end
        else
            lpCurr = lps(end);
        end
        
        
        %Set year in output:
        dateOut(ii,1) = dateRef(1) + nYrs;
        
        %Find month:
        if nPrec > 1
            if lpCurr
                dateOut(ii,2) = find(daysMnthCumLp > daysIn(ii), 1, 'first');
                if dateOut(ii,2) == 1
                    daysPrev = 0;
                else
                    daysPrev = daysMnthCumLp(dateOut(ii,2)-1);
                end
            else
                dateOut(ii,2) = find(  daysMnthCum > daysIn(ii), 1, 'first');
                if dateOut(ii,2) == 1
                    daysPrev = 0;
                else
                    daysPrev = daysMnthCum(dateOut(ii,2)-1);
                end
            end

            daysIn(ii) = daysIn(ii) - daysPrev;

            %Find day:
            if nPrec > 2
                dateOut(ii,3) = 1 + daysIn(ii);
            end
        end
    end 
end

%Find hours, minutes, seconds
if nPrec > 3
    for ii = 1 : numel(daysIn)
        for jj = 4 : nPrec
            switch jj
                case 4
                    divisor = 24;
                case 5 || 6
                    divisor = 60;
            end
            int = floor(dateOut(ii,jj-1));
            dateOut(ii,jj) = (dateOut(ii,jj-1) - int)*divisor;
            dateOut(ii,jj-1) = int;
        end
    end
end


%%Form Outputs using variable arguments
%If only one output argument and only one daysIn, then output is [year, month, day] 
if nargout <= 1
    varargout = {dateOut};
else
    nArg = nargout;
    
    if nArg ~= nPrec
       error('days2Date:nArgOut',['There are ' num2str(nArg) ' output arguments, but ' num2str(nPrec) ' are required.']);
    end
    
    varargout = cell(nArg,1);
    
    for jj = 1 : nArg
       varargout{jj} = dateOut(:,jj); 
    end
end





% for ii = 1 : numel(daysIn)
%     %Initialize cumulative days:
%     datCurr = 0;
%     
%     %If either 'daysIn' of the month/year of the reference date are NaN, output NaN:
%     if any([isnan(daysIn(ii)), isnan(dateRef(1:2)) ]) > 0
%         continue
%     end
%     
%     %Reset current ouput day to 1 (from daysRef) so that calculation 
%     % of days in year/month happens from the beginning each the year/month:
%     if numel(dateRef) >= 3
%         varargout{3}(ii) = 1;
%         daysIn(ii) = daysIn(ii) - dateRef(3) + 1;
%     end
%     
%     %Find output year:
%     while daysIn(ii) < daysIn(ii,:)
%         if blCal == 0 && is_leap_v2(varargout{1}(ii))
%             daysAddYr = 366;
%         elseif blCal == 1 || (blCal == 0 && ~is_leap_v2(varargout{1}(ii)))
%             daysAddYr = 365;
%         elseif blCal == 2
%             daysAddYr = 360;
%         else
%             error('days_2_date:unkownCal','Unkown calednar type input.');
%         end
%         
%         %Add days to sum and one year to output if cumulative sum less than
%         %or equal to daysIn
%         if daysIn(ii) + daysAddYr <=  daysIn(ii,:)
%         	daysIn(ii) = daysIn(ii) + daysAddYr;
%             varargout{1}(ii) = varargout{1}(ii) + 1;
%         else
%             break
%         end
%     end
% 
%     %Find month:
%     while daysIn(ii) < daysIn(ii,:)
% 
%         %May iterate past December, in which case must set to January 
%         if varargout{2}(ii) == 13
%            varargout{2}(ii) = 1;
%            varargout{1}(ii) = varargout{1}(ii) + 1;
%         end
%         
%         if varargout{2}(ii) == 2 && blCal == 1 %Select if February and 365 day calendar selected.
%             daysAddMnth = 28;
%         elseif blCal == 2
%             daysAddMnth = 30;
%         else
%             daysAddMnth = eomday(varargout{1}(ii),varargout{2}(ii));
%         end
% 
%         if daysIn(ii) + daysAddMnth <=  daysIn(ii,:)
%             daysIn(ii) = daysIn(ii) + daysAddMnth;
%             varargout{2}(ii) = varargout{2}(ii) + 1;
%         else
%             break
%         end
%     end
% 
%     
%     if varargout{2}(ii) > 12 %If month went over 12, transform to 1 
%         varargout{2}(ii) = varargout{2}(ii) - 12;
%         varargout{1}(ii) = varargout{1}(ii) + 1;
%     elseif varargout{2}(ii) == 0 %If month is 0, transform  1 (probably cannot occur)
%         varargout{1}(ii) = varargout{1}(ii) - 1;
%         varargout{2}(ii) = 12;
%     end
% 
% 
%     %Find day of month:
%     varargout{3}(ii) = daysIn(ii,:) - daysIn(ii) + 1;
% 
%     %Record day fraction if nargout = 4
%     if numel(varargout(:)) > 3
%         for jj = 4 : numel(varargout(:))
%             switch jj
%                 case 4
%                     divisor = 24;
%                 case 5 || 6
%                     divisor = 60;
%             end
%             int = floor(varargout{jj-1}(ii));
%             varargout{jj}(ii) = (varargout{jj-1}(ii) - int)*divisor;
%             varargout{jj-1}(ii) = int;
%         end
%     end
% end

