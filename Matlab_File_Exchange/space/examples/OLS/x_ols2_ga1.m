%This example script estimates some election data via  OLS and MESS-AR. 
%
%You can see the following for more information on MESS:
%
%LeSage, James and R. Kelley Pace, “Spatial Dependence in Data Mining,” 
%Data Mining for Scientific and Engineering Applications, 
%Edited by Robert L. Grossman, Chandrika Kamath, Philip Kegelmeyer, Vipin Kumar, and Raju R. Namburu, 
%Kluwer Academic Publishing, 2001.
%
%This uses the election data (dependent variable ln(voter turnout/voter population) examined in:
%
%Pace, R. Kelley, and Ronald Barry, Quick Computation of Regressions with a Spatially Autoregressive Dependent Variable,
%Geographical Analysis, Volume 29, Number 3, July 1997, p. 232-247. 
%
%Written by Kelley Pace, www.spatial-statistics.com, 12/25/02

clear all;
close all;
format short;
clc

%core independent variables
load xsub;
%dependent variables
load y;
%locational coordinates
load xcoord;
load ycoord;

var_labels=strvcat('Variables ' , 'Voting Pop ','Education ','Home Ownership ','Income ','Lag Voting Pop','Lag Education','Lag Home Ownership  ','Lag Income', 'Intercept','Alpha');


%number of observations
n=length(y);

%finding Delaunay weight matrices
[wswdel,wwsdel,wmatdel]=fdelw2(xcoord,ycoord);

%selecting a symmetric weight matrix
d=wswdel;

%forming matrix of independent variables
x=[xsub d*xsub ones(n,1)];

%number of observations, variables
[n,k]=size(x);


%Estimating OLS
tic
%[logliks_ols, emax_ols, bmax_ols, srds_ols, prhigher_ols]=fols2(x, y);
 [ bmax_ols, srds_ols, prhigher_ols, emax_ols, maxlik_ols]=fols2(x, y);
OLS_time=toc

%Estimating MESS-CAR
tic
%[alphamax, logliks, emax, bmax, srds, prhigher]=fmess_car2(x,y,d);
[bmax, srds, prhighers, emax, logliks]=fmess_car2(x, y, d);
MESS_time=toc

ests=[bmax_ols srds_ols bmax srds];

%using James LeSage's printing routine, www.spatial-econometrics.com
info.rnames=var_labels;
info.cnames=strvcat('OLS Beta ', 'OLS SRDs  ','MESS Beta ', 'MESS SRDs  ');
mprint(ests, info);
 
disp(' ')
disp('Obviously, the spatial parameter is quite significant given the signed root deviance (SRD) on alpha.');
disp('Note, the change in the income parameters and associated SRDs.');





 

 
   
 
 
 




