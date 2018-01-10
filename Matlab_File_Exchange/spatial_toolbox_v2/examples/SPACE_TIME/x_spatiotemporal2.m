clear all
clear global
close all
format bank
clc
%Written by Kelley Pace, www.spatial-statistics.com, 1/9/03.
%
%If you use the fmex temporal nearest neighbor routine, you must have tnmex4.dll in the search path.
%Naturally, this applies only to windows. The source code is available for other platforms.
%
%This example uses some Baton Rouge housing data as discussed in:
%
%Pace, R. Kelley and Ronald Barry, O.W. Gilley, C.F. Sirmans, 
%“A Method for Spatial-temporal Forecasting with an Application to Real Estate Prices,” 
% International Journal of Forecasting, Volume 16, Number 2, April-June 2000, p. 229-246.
%
%I have included these lines to document the variable definitions.
%
%y is mean-centered log(sales price) 
%
%xcoord, ycoord are rotated and translated latitude and longitude
%
%tdums are indicator variables based upon the date the house sold
%s1=saleyr==1985;
%s2=saleyr==1986;
%s3=saleyr==1987;
%s4=saleyr==1988;
%s5=saleyr==1989;
%s6=saleyr==1990;
%s7=saleyr==1991;
%s8=(saleyr==1992)+(saleyr==1993);
%tdums=[ s2 s3 s4 s5 s6 s7 s8 ];
%
%xquant=[log(age+1) log(livesf) log(oarea) (fullbath+halfbath)] where these variables are mean-centered
%age is the age of the house.
%livesf is the living area of the house
%oarea is other area in the house
%fullbath and halfbath refer to full and half bathrooms
%
%The observations are ordered by transaction date with the oldest ones in the first rows and the newest ones in the last rows.
%These transformations were done to preserve confidentiality.

load y;
load xquant;
load tdums;
load xcoord;
load ycoord;

%Some intial constant definitions

%number of intial observations to exclude in the spatial-temporal autoregression
firstobs=300;
%number of observations to use in computing temporal lags
mt=180;
%number of spatial neighbors
ms=30;
%geometric weighting factor given to neighbors
rho=0.80;
%first observation to use in computing spatiotemporal neighbors
istart=100;

%percentiles to be computed
pers=[0 1 5 10 25 50 75 90 95 99 100];

%number of observations
n=length(y);
%constant vector
o=ones(n,1);

%computation of previously transacted nearest neighbors
%tic;
%this function uses matlab only and works OK for moderate sample sizes under ver 6.5 on PCs.
%Earlier versions of Matlab may be very slow.
%nnmat1=fnn_ml2(xcoord, ycoord, ms, n, istart);
%mltime=toc

tic;
%This is a fortran 90 mex file and operates faster and can handle larger sample sizes.
%However, it still rises superlinearly with n.
nnmat=fnn_mex2(xcoord, ycoord, ms, n, istart);
%There are differences in the functions due to duplicate locations and lack of perfect precision.
temporal_nn_time=toc

%forming OLS independent variables
xolsall=[o ycoord xcoord xcoord.*ycoord xcoord.^2 ycoord.^2 tdums xquant];

%throwing away some initial observations for comparability with spatiotemporal method
xols=xolsall((firstobs+1):end,:);
yols=y((firstobs+1):end);

%Estimating via OLS
%[olslogliks, olsemax, olsbmax, olssrds, olsprhigher]=fols2(xols, yols);
 [olsbmax, olssrds, olsprhigher, olsemax, olslogliks]=fols2(xols, yols);
 
%need to strip out a trailing 0 used to match spatial routines
olsbmax=olsbmax(1:(end-1));
olssrds=olssrds(1:(end-1));

%Some OLS statistics
ydifols=yols-mean(yols);
ssqols=ydifols'*ydifols;
olssse=olsemax'*olsemax;
olsr2=1-olssse/ssqols;
olseper=prctile(olsemax,pers);
olsmaxlik=max(olslogliks);
[nols,kols]=size(xols);

%printing these out
disp('Spatial-Temporal Modeling of Housing Prices')
disp(' ')
disp('OLS Results with Time Dummies')
info1.rnames=strvcat('Variables', 'Intercept','lat','lon','lat*lon','lat*lat','lon*lon','I(86)','I(87)','I(88)','I(89)','I(90)','I(91)','I(92,93)','ln(age+1)','ln(Living Area)','ln(Other area)','Baths');
info1.cnames=strvcat('OLS b ','OLS SRDS');
info0.rnames=strvcat(' ','n','k','Loglik','R-squared');

mprint([olsbmax olssrds], info1)
mprint([nols kols olsmaxlik olsr2]',info0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spatiotemporal model

%index rising linearly with time
indext1=(0:(n-1))';
%east-west coordinate -- although the rotation has changed this meaning
ew=xcoord;
%north-south coordinate -- although the rotation has changed this meaning
ns=ycoord;
%a shorter variable name
z=xquant;
%number of observations
[n,k]=size(z);

%for dependent variable
%temporal lag of all ln market prices 
lagty=movav2(y,mt);
%spatial lag (although also previous in time) of ln market prices
lagdy=fts(y, ms, rho, nnmat);
%spatial lag of the temporal lag
lagdty=fts(lagty, ms, rho, nnmat);
%temporal lag of the spatial lag
lagtdy=movav2(lagdy,mt);

%same transformations for subset of independent variables
lagtz=movav2(z, mt);
lagdz=fts(z, ms, rho, nnmat);
lagdtz=fts(lagtz, ms, rho, nnmat);
lagtdz=movav2(lagdz, mt);

%A general specification involving all types of lags
%xallt=[o indext1 ns ew  (z-lagtz) (lagdz-lagdtz) lagtz lagdtz  lagtdz  lagdy lagty  lagdty lagtdy];

%A subset of the general specification
xallt=[o    (z-lagtz) (lagdz-lagdtz)  lagdy   lagdty  lagtdy];

%discarding initial observations (like dropping initial lag in time series)
xst=xallt( ((firstobs+1):end) ,:);
%time differencing dependent variable
yallt=y-lagty;
%discarding initial observations (like dropping initial lag in time series)
yst=yallt((firstobs+1):end);

%Applying OLS -- can do this since all conditioning is on past values -- consistent
%[logliks, emax, bmax, srds, prhigher]=fols2(xst, yst);
 [bmax, srds, prhigher, emax, logliks]=fols2(xst, yst);
 bmax=bmax(1:(end-1));%need to strip out a trailing 0 used to match spatial routines
 srds=srds(1:(end-1));%need to strip out a trailing 0 used to match spatial routines
 

%usual statistics
sse=emax'*emax;
ystdif=yst-mean(yst);
ssq=ystdif'*ystdif;
r2=1-sse/ssq;
eper=prctile(emax,pers);
maxlik=max(logliks);
[nspace,kspace]=size(xst);

 %printing results
 disp(' ')
 disp('Dependent Variable: (I-T)y')
 info2.rnames=strvcat('Variables', 'Intercept', '(I-T)ln(age+1)','(I-T)ln(Living Area)','(I-T)ln(Other area)','(I-T)Baths','D(I-T)ln(age+1)','D(I-T)ln(Living Area)','D(I-T)ln(Other area)','D(I-T)Baths','Dy','DTy','TDy');
info2.cnames=strvcat('Space-time b ','Space-time SRD');
mprint([bmax srds], info2)
info3.rnames=strvcat(' ','mt','rho','n','k','Loglik','R-squared');
mprint([mt rho nspace kspace maxlik r2]',info3)
disp(' ')
disp('Note, the improvement in the likelihood using fewer variables than the trend surface, time dummy model.')
disp(' ')

%computing recursive estimates (starting with observation 1000)
[erecur, sset, brecur]=frecursive2(xst,yst,1000);

%graphing 1-out errors
disp(' ')
disp('The recursive residuals, which are from 1-out predictions show nice random behavior')
plot(erecur(1000:end),'.')
xlabel('Index')
ylabel('Recursive Residuals')
title('Recursive Residuals over Time')

