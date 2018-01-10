%This uses the election data (dependent variable ln(voter turnout/voter population) examined in:
%
%Pace, R. Kelley, and Ronald Barry, Quick Computation of Regressions with a Spatially Autoregressive Dependent Variable,
%Geographical Analysis, Volume 29, Number 3, July 1997, p. 232-247. 
%
%Written by Kelley Pace, www.spatial-statistics.com, 12/24/02

clear all;
clear global;
clc;
close all;
%core independent variables
load xsub;
%dependent variables
load y;
%locational coordinates
load xcoord;
load ycoord;

info.cnames=strvcat('Beta Estimates ', 'Signed Root Deviances  ', 'PR of Higher SRDS ');
var_names=strvcat('Variables ' , 'Voting Pop ','Education ','Home Ownership ','Income ', 'Intercept','Alpha');
info.rnames=var_names;

%finding Delaunay weight matrices
tic;
[wswdel,wwsdel,wmatdel]=fdelw2(xcoord,ycoord);
weight_matrix_time=toc;

%selecting a symmetric weight matrix
d=wswdel;
sym=1;

%number of observations
n=length(y);

%forming matrix of independent variables
x=[xsub  ones(n,1)];

%Estimating a mixed regressive spatially autoregressive model without log-determinant 
%using Chebyshev approximations, Taylor bounds, and likelihood dominance inference
tic
 [bmax1, srds1, prhighers1, emax1, logliks1]=fcar2(x, y, d);
time_sans_lndet=toc;

ests1=[bmax1 srds1 prhighers1];

disp(' ')
disp('fcar2 Estimation Results using Chebyshev lndet approximation and likelihood dominance inference');
mprint(ests1, info)
disp(' ')


%computing log-determinants
tic;
[detvalz]=fdet_interp2(d, sym);
log_det_time=toc;
 


%number of observations, variables
[n,k]=size(x);

%Estimating CAR
tic
 [bmax, srds, prhighers, emax, logliks]=fcar2(x, y, d, detvalz);
time_with_lndet=toc;

ests=[bmax srds prhighers];

disp('fcar2 Estimation Results using exact lndet');
mprint(ests, info)
disp(' ')
disp('Note, the likelihood dominance SRDs are smaller in magnitude than the exact SRDs.');
disp('However, they can still document statistical significance for many variables, and thus can');
disp('prove useful in many circumstances.');
disp(' ')
disp(' ')


tic
 [bmax2, srds2, prhighers2, emax2, logliks2]=fcar2(x, y, d, detvalz,100);
time_with_lndet_addlpoints=toc;

ests2=[bmax srds prhighers];
disp('fcar2 Estimation Results using exact lndet and additional interpolation points');
mprint(ests2, info)
disp(' ')

time_vec=[weight_matrix_time;log_det_time;time_sans_lndet;time_with_lndet;time_with_lndet_addlpoints];
info1.rnames=strvcat('Times','Weight Matrix','Lndet','Approximate CAR','Exact CAR', 'Exact CAR +');
disp('Times (must add lndet time to exact CAR or exact CAR with added points to compare to approximate CAR)');
mprint(time_vec,info1)
disp(' ')
disp('The cost of going from 20 to 100 additional beta computations is not material for this example.');
disp('The advantages of interpolation become clearer for large k. In this example, the additional');
disp('interpolation points do not add more accuracy, because the CAR dependence parameter is at its');
disp('maximum value. This is not abnormal behavior for CAR.');
disp(' ')
disp(' ')

plot(detvalz(:,1),logliks)
title('Profile likelihoods vs \alpha for global model and delete-1 submodels');
xlabel('Dependence parameter (\alpha)');
ylabel('Profile log-likelihood');
legend(strvcat('Global likelihood',var_names(2:end-1,:)),-1);
disp('This function also produces profile likelihoods across a range of values of the dependence parameter,');
disp('as shown in the figure.');




 






 

 
   
 
 
 




