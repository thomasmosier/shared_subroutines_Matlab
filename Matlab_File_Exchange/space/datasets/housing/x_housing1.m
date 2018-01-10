%This example program uses the data from:
%
%Pace, R. Kelley, and James LeSage, “Chebyshev Approximation of Log-determinants of Spatial Weight Matrices,” 
%Computational Statistics and Data Analysis, forthcoming.
%
%
%The functions used here are described in
%
%Pace, R. Kelley, and Ronald Barry, "Quick Computation of Regressions with a Spatially Autoregressive Dependent Variable,"
%Geographical Analysis, Volume 29, Number 3, July 1997, p. 232-247. 
%
%Pace, R. Kelley, and Dongya Zou, 
%“Closed-Form Maximum Likelihood Estimates of Nearest Neighbor Spatial Dependence,” 
%Geographical Analysis, Volume 32, Number 2, April 2000, p. 154-172.
%
%Pace, R. Kelley, and James LeSage, “Likelihood Dominance Spatial Inference,” forthcoming Geographical Analysis.
%
%Pace, R. Kelley, and James LeSage, “Chebyshev Approximation of Log-determinants of Spatial Weight Matrices,” 
%Computational Statistics and Data Analysis, forthcoming.
%
%LeSage, James and R. Kelley Pace, “Spatial Probit and Tobit,” Spatial Statistics and Spatial Econometrics, Edited by Art Getis, Palgrave, 2003.
%
%Written by Kelley Pace, www.spatial-statistics.com,  8/98, revised %12/26/02.
%
%Written by Kelley Pace, www.spatial-statistics.com, 1/11/03.

clear all;
clear global;
close all;
format short;
%clc

%locations on a plane
load xcoord;
load ycoord;
load xsub;
load y;

[n, psub]=size(xsub);
o=ones(n,1);

%finds the very closest neighbor - supplied as a list
tic;
 [nnlist]=fclosest_nn2(xcoord, ycoord);
closest_neighbor_search_time=toc;

%finding Delaunay weight matrices
tic;
[wswdel,wwsdel,wmatdel]=fdelw2(xcoord,ycoord);
delaunay_time=toc;

x=[xsub xsub(nnlist,:) o ];

%%%%%%%%%%%%%%%%%%%%% Closest neighbor %%%%%%%%%%%%%%%%%%%%%%%

%computes mixed estimates using the list of nearest neighbor indices
tic;
[cbmax, csrds, cprhigher,cemax, cmaxlik]=fclosest_ar2(x, y, nnlist);
closest_ar_time=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% Comparing to MESS %%%%%%%%%%%%%%%%%%%%%%

%selecting a symmetric weight matrix
tic;
d=wswdel;%if you don't wish to have a doubly stochastic weight matrix
%d=fdoubly2(wswdel);%if you wish to have a doubly stochastic weight matrix
doubly_stochastic_time=toc;

%forming matrix of independent variables
x=[xsub d*xsub o];

%Estimating MESS-AR
tic
[bmax3, srds3, prhigher3, emax3, maxlik3]=fmess_ar2(x, y, d);
MESS_time=toc;

delaunay_loglik=max(maxlik3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Comparing to Mix with Likelihood Dominance %%%%%%%%%%%%

%Estimating a mixed regressive spatially autoregressive model without log-determinant 
%using Chebyshev approximations, Taylor bounds, and likelihood dominance inference
tic
 [bmax4, srds4, prhighers4, emax4, logliks4]=fmix2(x, y, d);
mix_time_sans_lndet=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Comparing to OLS with and without lagged independent variables %%%%%%%%%%%%

%completely aspatial OLS
tic;
[bmax_aspatial, srds_aspatial, prhigher_aspatial, emax_aspatial, maxlik_aspatial]=fols2([xsub o], y);
aspatial_ols_time=toc;

%putting this in a vector to match the lengths of the other routines
bols_aspatial=zeros(length(bmax4),1);
bols_aspatial(1:psub)=bmax_aspatial(1:psub);
bols_aspatial(end-1)=bmax_aspatial(psub+1);

tols_aspatial=zeros(length(srds4),1);
tols_aspatial(1:psub)=srds_aspatial(1:psub);
tols_aspatial(end-1)=srds_aspatial(psub+1);

aspatial_loglik=max(maxlik_aspatial);

%OLS with spatial lags of independent variables
tic;
[ bmax_ols, srds_ols, prhigher_ols, emax_ols, maxlik_ols]=fols2(x, y);
ols_time=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Printing Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info.width=132;
info.rnames=strvcat('Variables ' , 'Land area ','Pop ','Per cap Income ','Age ','Lag Land area','Lag Pop','Per cap Income  ','Lag Age', 'Intercept','Alpha');
info.cnames=strvcat('b OLS   ', 'b closest    ', 'b MESS ', 'b Mix ');
info1.rnames=info.rnames;
info1.cnames=strvcat('SRDs OLS ', 'SRDs closest ', 'SRDs MESS ', 'SRDs Mix ');

disp('Comparing OLS, closest neighbor (NN), MESS, Mixed in Beta, dependence estimates');
disp(' ');
mprint([bols_aspatial  cbmax  bmax3  bmax4],info)
disp(' ');
disp('Comparing OLS, closest neighbor (NN), MESS, approximate Mixed in SRDS');
disp(' ');
mprint([ tols_aspatial  csrds  srds3 srds4],info1)
disp(' ');
disp('The closest neighbor approach is intermediate to a non-spatial approach (OLS) and a full spatial approach (MESS)');
disp('or the approximate mixed routine. Note, the close agreement between MESS and the approximate mixed routine.');
disp('Keep in mind that the approximate mixed routine uses likelihood dominance inference which is a lower bound to');
disp('the signed root deviance (SRD). Also, OLS in this case uses the spatial averages of the basic independent variables as');
disp('additional independent variables.');
disp(' ');


tmat=[ols_time;closest_ar_time;MESS_time;mix_time_sans_lndet;doubly_stochastic_time];
disp(' ' );
disp('The following table shows the timings in seconds for the routines -- all are very fast.');
disp(' ');

%making a table using James LeSage's mprint routine at
%www.spatial-econometrics.com
info2.rnames=strvcat('ESTIMATION TIME TABLE: ', '  OLS','  Closest AR', '  MESS', '  Approximate Mix','  Doubly Stochastic Scaling');
mprint(tmat,info2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% NN computations %%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(' ')
disp('Just selecting a particular weight matrix seems arbitrary. Here we take 30 nearest neighbors');
disp('and weight these geometrically. A rho of 1 indicates no decline in the weight given to further');
disp('neighbors relative to closer ones, while a rho of 0.5 would give half the weight to the 2nd NN');
disp('as it would to the 1st NN. Thus, rho allows changes in the effective number of neighbors used');
disp('without actually varying the number of neighbors. It often makes sense in this approach to set');
disp('the number of neighbors to a fairly high level (such as 30)');
disp(' ');

%needs delaunay weights
%Computes nearest neighbors and records timings
nnmax=30;
tic;
 fneighbors2(wswdel, xcoord, ycoord, nnmax);
neighbor_time=toc;

%Setting the number of neighbors to the maximum and iterating over a grid of geometric
%weighting parameters to find largest likelihood

%weighting parameter grid
rhovec=[0.8 0.85 0.9 0.95 1];
%number of elements in the grid
rhoiter=length(rhovec);

%a loop to find the maximum likelihood across rhovec using MESS
tic;
rez=zeros(rhoiter,1);
for i=1:rhoiter
rho=rhovec(i);
%Symmetricizes and weights NNs 
[wswnn,wwsnn,wmatnn]=fsym_neighbors2(nnmax, rho);
%d=fdoubly2(wswnn);%if you wish for doubly stochastic scaling;
d=wswnn;%if you don't want doubly stochastic scaling
[bmaxi, srdsi, prhigheri, emaxi, maxliki]=fmess_ar2(x, y, d);
rez(i)=max(max(maxliki));
betarez(:,i)=bmaxi;
srdsrez(:,i)=srdsi;
end;
rho_search_time=toc;

rezmat=[rhovec' rez];
info3.cnames=strvcat('rho', '   likelihood');
mprint(rezmat,info3);

[maxlikrho,maxindrho]=max(rez);
optrezvec_spatial=rezmat(maxindrho,:)';
optrezvec=[optrezvec_spatial; delaunay_loglik;aspatial_loglik];
disp(' ');
info3a.rnames=strvcat(' ','  optimum rho', '  maximum likelihood across rho','  Delaunay maximum likelihood','  aspatial likelihood (OLS)');
mprint(optrezvec,info3a);
disp(' ');
disp('Examining the MESS loglikelihoods over rho and Delaunay and contrasting it with the loglikelihood from');
disp('appling OLS to the basic non-spatial independent variables shows that even a suboptimal');
disp('choice of rho or Delaunay still dominates the use of an aspatial model in this case and that optimizing over rho dominates');
disp('an arbitrary choice of weight matrices.');
disp(' ');
disp(' ');

disp('It did not take long to optimize over some possible weight matrices');
disp(' ' );
disp('The following table shows the timings for computing the nearest neighbors and optimizing over rho')
disp(' ');
%making a table using James LeSage's mprint routine at
%www.spatial-econometrics.com

info4.rnames=strvcat('TIME TABLE: ', '  NN computation','  Time to find optimum rho');
tmat2=[neighbor_time;rho_search_time];
mprint(tmat2,info4);

disp(' ');
disp(' ');
disp('MESS results for optimum rho');
disp(' ');
info5.rnames=info.rnames;
info5.cnames=strvcat('Optimal b MESS', 'Optimal SRDS MESS');
bmaxopt=betarez(:,maxindrho);
srdsopt=srdsrez(:,maxindrho);
optmess=[bmaxopt srdsopt];
mprint(optmess,info5);

