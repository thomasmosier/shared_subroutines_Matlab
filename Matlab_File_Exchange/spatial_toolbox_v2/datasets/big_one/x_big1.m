%This example script simulates a MESS-AR process, and estimates it. 
%
%You can see the following for more information:
%LeSage, James and R. Kelley Pace, “Spatial Dependence in Data Mining,” 
%Data Mining for Scientific and Engineering Applications, 
%Edited by Robert L. Grossman, Chandrika Kamath, Philip Kegelmeyer, Vipin Kumar, and Raju R. Namburu, 
%Kluwer Academic Publishing, 2001.
%
%Written by Kelley Pace, www.spatial-statistics.com, 12/24/02
%
%Don't use 1M unless you have some memory in your machine (more than 512mb and probably 1gb)!
%
clear global;
clear all;
close all;
clc

%number of observations
n=1000000;

%number of non-constant independent variables
p=4;

%number of Monte Carlo iterations
iter=1;

%Standard deviation of the error term
sigma=1.00;

%Value of the spatial dependence parameter in the matrix exponential
%covariance function. If this is high, say 6 or over, one needs to
%increase q in the simulation and estimation routines.
alpha=4;

%generating artificial locational coordinates
xcoord=rand(n,1);
ycoord=rand(n,1);

%finding Delaunay weight matrices
tic;
[wswdel,wwsdel,wmatdel]=fdelw2(xcoord,ycoord);
delauany_time=toc

%selecting a symmetric weight matrix
d=wswdel;
clear wswdel wwsdel wmatdel

%creating random independent variables and including an intercept
x=[randn(n,p) ones(n,1)];
[n,k]=size(x);

%setting the true regression parameters
beta=ones(k,1);

%true predictions and the random error term
rv=[x*beta randn(n, iter)];

%simulating exp(alpha/2)*X*beta, exp(alpha/2)*error
tic;
rvs=fmess_sim2(d,rv,alpha);
simulation_time=toc

%the simulated dependent variable is exp(alpha/2)*X*beta+sigma*exp(alpha/2)*error 
ymat=rvs(:,1)+rvs(:, 2) *sigma;

%Estimating MESS-AR
tic
    [bmaxs, srds, prhigher, emax, maxlik]=fmess_ar2(x, ymat, d);
estimation_time=toc
disp('finished estimation')


 
 %true and estimated regression and spatial dependence parameters
disp(' ')
disp(' ')
 %Using James LeSage's printing function, www.spatial-econometrics.com
 info.cnames=strvcat('true beta ','mean beta MESS ');
 mprint([[beta;alpha] bmaxs ],info)
 
 disp(' ')
 disp('The truth and estimates agree nicely for large samples.')
 

 
  
 
 





 

 
   
 
 
 




