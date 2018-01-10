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

info.cnames=strvcat('Beta Estimates ', 'Signed Root Deviances  ', 'PR of Higher SRDS ');
var_names=strvcat('Variables ' , 'Voting Pop ','Education ','Home Ownership ','Income ','Lag Voting Pop','Lag Education','Lag Home Ownership  ','Lag Income', 'Intercept','Alpha');
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
x=[xsub d*xsub ones(n,1)];

%Estimating a mixed regressive spatially autoregressive model without log-determinant 
%using Chebyshev approximations, Taylor bounds, and likelihood dominance inference
tic
 [bmax1, srds1, prhighers1, emax1, logliks1]=fmix2(x, y, d);
mix_time_sans_lndet=toc;

ests1=[bmax1 srds1 prhighers1];

disp('fmix2 Estimation Results using Chebyshev lndet approximation and likelihood dominance inference');
mprint(ests1, info)
disp(' ')


%computing log-determinants
tic;
[detvalz]=fdet_interp2(d, sym);
log_det_time=toc;
 


%number of observations, variables
[n,k]=size(x);

%Estimating a mixed regressive spatially autoregressive model
tic
 [bmax, srds, prhighers, emax, logliks]=fmix2(x, y, d, detvalz);
mix_time_with_lndet=toc;

ests=[bmax srds prhighers];

disp('fmix2 Estimation Results using exact lndet');
mprint(ests, info)
disp(' ')
disp('Note, the likelihood dominance SRDs are smaller in magnitude than the exact SRDs.');
disp('However, they can still document statistical significance for many variables, and thus can');
disp('prove useful in exploration.');
disp(' ')

time_vec=[weight_matrix_time;log_det_time;mix_time_sans_lndet;mix_time_with_lndet];
info1.rnames=strvcat('Times','Weight Matrix','Lndet','Approximate Mix','Exact Mix');
disp('Times (must add lndet time to exact mix to compare to approximate mix)');
mprint(time_vec,info1)
disp(' ')

plot(detvalz(:,1),logliks)
title('Profile likelihoods vs \alpha for global model and delete-1 submodels');
xlabel('Dependence parameter (\alpha)');
ylabel('Profile log-likelihood');
legend(strvcat('global lik',var_names(2:end-1,:)),-1);
disp('This function produces profile likelihoods across a range of values of the dependence parameter.');




 






 

 
   
 
 
 




