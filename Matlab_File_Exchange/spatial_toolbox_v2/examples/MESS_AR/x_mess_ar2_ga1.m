%This example script simulates a MESS-AR process, and estimates it. 
%
%You can see the following for more information:
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
%Written by Kelley Pace, www.spatial-statistics.com, 12/24/02

clear all;
close all;
clc;

%core independent variables
load xsub;
%dependent variables
load y;
%locational coordinates
load xcoord;
load ycoord;

info.rnames=strvcat('Variables ' , 'Voting Pop ','Education ','Home Ownership ','Income ','Lag Voting Pop','Lag Education','Lag Home Ownership  ','Lag Income', 'Intercept','Alpha');
info.cnames=strvcat('Beta Estimates ', 'Signed Root Deviances  ', 'PR of Higher SRDS ');

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

%Estimating MESS-AR
tic
[bmax, srds, prhigher, emax, maxlik]=fmess_ar2(x, y, d);
MESS_time=toc

ests=[bmax srds prhigher];
mprint(ests, info)
 




 

 
   
 
 
 




