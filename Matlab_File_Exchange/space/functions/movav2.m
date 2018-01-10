function mavz=movav2(z,mt);
%
%mavz=movav2(z,mt);
%
%This function computes a moving average of the past mt observations of z.
%The function uses less than mt observations for the first mt observations
%of z.
%
%INPUT:
%
%The n by k matrix z
%
%The scalar mt which gives the past window
%
%
%OUTPUT:
%
%mavz is a n by k matrix of the smoothed z
%
%NOTES:
%
%This function was used in:
%
%Pace, R. Kelley and Ronald Barry, O.W. Gilley, C.F. Sirmans, 
%“A Method for Spatial-temporal Forecasting with an Application to Real Estate Prices,” 
%International Journal of Forecasting, Volume 16, Number 2, April-June 2000, p. 229-246.
%
%Written by Kelley Pace, www.spatial-statistics.com, 12/24/02

[n, k]=size(z);
cumz=cumsum(z);
qseqinv=1./(1:mt)';
scalemat=qseqinv(:,ones(1,k));
ccc1=cumz(1:mt,:).*scalemat;
ccc2=(cumz((mt+1):(n-1),:)-cumz(1:(n-mt-1),:))/mt;
mavz=[z(1,:);ccc1;ccc2];


