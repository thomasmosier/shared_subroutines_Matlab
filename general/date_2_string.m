function strDate = date_2_string(date)

if ~any(size(date) == 1)
   error('date2str:multipledates', 'Multiple dates were passed to this function, which has not been programmed for.') 
end

strDate = '';
nPrec = numel(date);
for ii = 1 : nPrec
    strDate = [strDate num2str(date(ii))];
    if ii ~= nPrec
        strDate = [strDate '-'];
    end
end