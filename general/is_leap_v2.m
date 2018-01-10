function blLeap = is_leap_v2(yr)

blLeap = zeros(size(yr), 'single');

for ii = 1 : numel(yr)
    if mod(yr(ii),4) == 0 && mod(yr(ii),100) ~= 0
        blLeap(ii) = 1;
    elseif mod(yr(ii),100) == 0 && mod(yr(ii),400) == 0
        blLeap(ii) = 1;
    end
end

%for ii = indChk

