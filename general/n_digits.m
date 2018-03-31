function y = n_digits(x)

x = abs(x); %in case of negative numbers
y = 0;

xp = x;

cntr = 0;
while (floor(xp)~= round2(x,y+1))
    y = y+1;
    xp = x*10;
    
    cntr = cntr + 1;
    
    if cntr > 10^6
       break 
    end
end

%https://www.mathworks.com/matlabcentral/answers/52068-machine-precision-problem-in-code

