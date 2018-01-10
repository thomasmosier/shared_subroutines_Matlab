function y = n_digits(x)

x = abs(x); %in case of negative numbers
y = 0;

xp = x;

while (floor(xp)~= round2(x,y+1))
    y = y+1;
    xp = x*10;
end

%https://www.mathworks.com/matlabcentral/answers/52068-machine-precision-problem-in-code

