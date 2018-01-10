%This script demonstrates Spatial Autoregressive Local Estimation as
%described in:
%
%Pace, R. Kelley, and James LeSage, “Spatial Autoregressive Local Estimation,” 
%Spatial Statistics and Spatial Econometrics, Edited by Art Getis, Palgrave, 2003.
%
%Using election data on 3,107 counties as described in:
%
%Pace, R. Kelley, and Ronald Barry, Quick Computation of Regressions with a Spatially Autoregressive Dependent Variable,
%Geographical Analysis, Volume 29, Number 3, July 1997, p. 232-247. 
%
%
%Written by Kelley Pace, www.spatial-statistics.com, 12/21/02
%
%This takes 3.4 minutes to run on a AMD 1700+ machine when using the 3,107
%observation election data for nsub=250. Takes 6.33 minutes for nsub=150 for 1 ghz PIII. 
%
clear all;
close all;
format short;
clc

%load dependent variable
load y;
%load non-constant independent variables
load xsub;
%loading geographical coordinates
load xcoord;
load ycoord;

%number of observations
n=length(y);
%constant term
o=ones(n,1);

%maximum number of observations in the local regressions
%I have set this to 150 for demonstration purposes -- lower settings run
%faster, and higher settings are better for examining the
%stability-prediction trade-offs. The variable nsub must be greater at
%least 50  due to the smoother routines used in some of the graphs below.
nsub=150;

%global delaunay weight matrix
[wswdel,wwsdel,wmatdel]=fdelw2(xcoord,ycoord);
d=wswdel;

%selecting expansion points
obsindices=1:3107;

%independent variables
xall=[xsub wwsdel*xsub o];
yall=[y wwsdel*y];

%number of observations, variaables
[n,k]=size(xall);

%minimum number of observations in the local regressions
start=2*k;

tic
%SALE function for multiple locii
[eouts, erecs, alphamaxes, bmaxes]=fwholesale2(xall, yall, d,  start, nsub, obsindices, xcoord, ycoord);
disp('finished SALE estimation')
toc
disp(' ')

%examining sequence of recursive residuals
es=median(abs(erecs(start:end,:)'))';
avmedy=movav2(es,25);
plot(start:nsub,avmedy,'.')
title('Smoothed SALE Recursive Residuals of Fringe Observations');
xlabel('Number of Local Observations');
ylabel('Absolute Error')
eval(['print -depsc -r1200 ' 'fig' num2str(1)]);


figure

%examining sequence of initial observation residuals
es2=median(abs(eouts(start:end,:)'))';
avmedy2=movav2(es2,25);
plot(start:nsub,avmedy2,'.')
title('Smoothed SALE Initial Holdout Residuals');
xlabel('Number of Local Observations');
ylabel('Absolute Error')
eval(['print -depsc -r1200  ' 'fig' num2str(2)]);

figure
%examining spatial dependence parameters across sample sizes
plot(median(alphamaxes((start+1):end,:)'), '.')
title('SALE Autoregressive Parameter Estimates')
xlabel('Number of Local Observations');
ylabel('Median Autoregressive Parameter Estimates')
eval(['print -depsc -r1200 ' 'fig' num2str(3)]);

ptiles=[0 10 25 50 75  90 100];

figure
a=alphamaxes(end,:)';
fdotplot2(xcoord, ycoord, a)
xlabel('Longitude')
ylabel('Latitude')
title('\alpha')
eval(['print -depsc -r1200 ' 'fig' num2str(4)]);

bdist=zeros(k,length(ptiles));
for ii=1:k
b=squeeze(bmaxes(end, ii, :))';
%bdist(ii,:)=[prctile(b, ptiles)]; If you have mathworks stat toolbox
bdist(ii,:)=[perctile(b, ptiles)'];%Peter Acklam's function pertile
figure
fdotplot2(xcoord, ycoord, b)
xlabel('Longitude')
ylabel('Latitude')
title(['\beta_{' num2str(ii) '}'])
%eval(['print -depsc -r1200 ' 'fig' num2str(4+ii)]);%optional saving of the graphics
end

%['Order Statistics for alpha, beta']
%orderstats=[ptiles;[prctile(a, ptiles)];bdist]; If you have mathworks stat toolbox
orderstats=[ptiles;[perctile(a, ptiles)'];bdist];%Peter Acklam's function pertile

%Printing using LeSage's mprint function (www.spatial-econometrics.com)
info.width=256;
info.rnames=strvcat(' ','Percentiles ' ,'alpha ', 'Voting Pop ','Education ','Home Ownership ','Income ','Lag Voting Pop','Lag Education','Lag Home Ownership  ','Lag Income', 'Intercept');
mprint(orderstats,info)

