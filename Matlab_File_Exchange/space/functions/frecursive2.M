function [erecur,sset,brecur]=frecursive2(xt,yt,start)
%
%[erecur,sset,brecur]=frecursive2(xt,yt,start)
%
%This function estimates recursive residuals, sum-of-squared errors,
%regression parameters.
%
%INPUT:
%
%xt represents n observations on k independent variables
%
%yt represents n observations on the dependent variable and its spatial
%lag. 
%
%start gives the minimum number of observations to use before beginning
%recursion
%
%OUTPUT:
%
%erecur represents nsub recursive residuals from the prediction of yt.
%Hence, erecur is a n by 1 vector.
%
%sset is a sequence of sum-of-squared residuals from y.
%
%brecur is the n by k recursive regression parameter estimates from
%estimating using y.
%
%NOTES:
%
%The values of erecur, sset, and brecur are zeros before observation start. The variable start must be greater than
%the number of variables. Recursive least squares is not the most numerically robust algorithm in the world, so be careful
%when using it with many correlated independent variables.
%
%This function was used in:
%
%Pace, R. Kelley and Ronald Barry, O.W. Gilley, C.F. Sirmans, 
%“A Method for Spatial-temporal Forecasting with an Application to Real Estate Prices,” 
% International Journal of Forecasting, Volume 16, Number 2, April-June 2000, p. 229-246.
%
%Written by Kelley Pace, www.spatial-statistics.com, 12/21/02.

%This begins by computing least squares for the first start-1 observations.
[n,k]=size(xt);
xtsub=xt(1:(start-1),:);
ytsub=yt(1:(start-1),:);
momtsub=xtsub'*xtsub;
invmomtsub=inv(momtsub);
btsub=invmomtsub*(xtsub'*ytsub);
etsub=ytsub-xtsub*btsub;

%initializing variables to be used in the recursive algorithm.
bj=btsub;
sj=invmomtsub;
momt=momtsub;
y2t=ytsub'*ytsub;

brecur=zeros(n,k);
dj=zeros(n,1);
erecur=zeros(n,1);
sset=zeros(n,1);

%This loop is the standard recursive least squares algorithm.
for i=start:n;

bj1=bj;
sj1=sj;

xj=xt(i,:)';
yj=yt(i,:)';
prej=xj'*bj1;
ej=yj'-prej;
erecur(i,:)=ej;

momt=momt+xj*xj';
y2t=y2t+yj*yj';

part1=sj1*xj;
part2=(1+part1'*xj);

dj(i)=sqrt(part2);
bj=bj1+(part1*ej)/part2;
sj=sj1-(part1*part1')/part2;
brecur(i,:)=bj';
sset(i)=y2t-bj'*momt*bj;
end;



