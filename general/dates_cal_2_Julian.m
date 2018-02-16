function [dayJ, datesOut] = dates_cal_2_Julian(datesIn, tRes)

%Calculate Julian dates:
if regexpbl(tRes,{'day','daily'}) %Calculate for day:
    if numel(datesIn(1,:)) == 3
        dayJ = Julian_day(datesIn,'ofyear');
        datesOut = datesIn;
    else
        error('toa_rad_DeWalle:dateNumel',['The date vector has ' num2str(numel(datesIn(1,:))) ' elements but 3 are needed.']);
    end
elseif regexpbl(tRes,'month') %if month time-step, calculate PET for each day and average to monthly value
    dayJ = nan(31*numel(datesIn(:,1)));
    datesOut = nan([31*numel(datesIn(:,1)),3]);
    cntr = 0;
    for ii = 1 : numel(datesIn(:,1))
        for jj = 1 : eomday(datesIn(ii,1),datesIn(ii,2))
            cntr = cntr + 1;
            datesOut(cntr,:) = [datesIn(ii,1:2),jj];
            dayJ(cntr) = Julian_day(datesOut(cntr,:),'ofyear');
        end
    end
    dayJ = dayJ(1:cntr);
    datesOut = datesOut(1:cntr,:);
else
    error('toa_rad_DeWalle:tUnit',['Current time units are ' tRes ' but days or months are required.  Program for this.'])
end