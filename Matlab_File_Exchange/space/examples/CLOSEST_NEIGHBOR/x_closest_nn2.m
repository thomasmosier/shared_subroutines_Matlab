%This example program uses the data from:
%
%Pace, R. Kelley, and Ronald Barry, "Quick Computation of Regressions with a Spatially Autoregressive Dependent Variable,"
%Geographical Analysis, Volume 29, Number 3, July 1997, p. 232-247. 

%The functions used here are described in
%
%Pace, R. Kelley, and Dongya Zou, 
%“Closed-Form Maximum Likelihood Estimates of Nearest Neighbor Spatial Dependence,” 
%Geographical Analysis, Volume 32, Number 2, April 2000, p. 154-172.
%
%Written by Kelley Pace, www.spatial-statistics.com,  8/98, revised %12/26/02.
%
%The function fclosestnn1 finds the very nearest (closest) neighbor to each observation
%The function fclosestmix1 computes the maximum liklelihood estimates given spatial dependence on the very nearest neighbor.

%Written by Kelley Pace, 8/98

clear all;
clear global;
close all;
format short;
clc

%locations on a plane
load xcoord;
load ycoord;
load xsub;
load y;

n=length(y);
o=ones(n,1);

%finds the very closest neighbor - supplied as a list
tic;
 [nnlist]=fclosest_nn2(xcoord, ycoord);
closest_neighbor_search_time=toc

%finding Delaunay weight matrices
tic;
[wswdel,wwsdel,wmatdel]=fdelw2(xcoord,ycoord);
delaunay_time=toc

x=[xsub xsub(nnlist,:) o ];

%computes mixed estimates using the list of nearest neighbor indices
tic;
[bmax, srds, prhigher,emax, maxlik]=fclosest_ar2(x, y, nnlist);
closest_ar_time=toc;


%loads nearest neighbor matrix -- the closest
load s1;

x=[xsub s1*xsub o];

%computes mixed estimates using the weight matrix s1
[bmax2, srds2, prhigher2emax2, maxlik2]=fclosest_ar2(x, y, s1);

disp('Note, the autoregressive parameter agree from the two different ways of specifying the closest neighbor.');
[bmax(end) bmax2(end)]

%%%%%%%%%%%%%% Comparing to MESS %%%%%%%%%%%%%%%%%%%%%%

%selecting a symmetric weight matrix
d=wswdel;

%forming matrix of independent variables
x=[xsub d*xsub o];

%Estimating MESS-AR
tic
[bmax3, srds3, prhigher3, emax3, maxlik3]=fmess_ar2(x, y, d);
MESS_time=toc;

%%%%%%%%%%%%%% Comparing to OLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
[ bmax_ols, srds_ols, prhigher_ols, emax_ols, maxlik_ols]=fols2(x, y);
ols_time=toc;

info.width=132;
info.rnames=strvcat('Variables ' , 'Voting Pop ','Education ','Home Ownership ','Income ','Lag Voting Pop','Lag Education','Lag Home Ownership  ','Lag Income', 'Intercept','Alpha');
info.cnames=strvcat('b OLS ', 'SRD OLS ', 'b NN ', 'SRD NN ','b MESS ', 'SRD MESS ');

disp('Comparing OLS, closest neighbor (NN), and MESS')
disp(' ')
mprint([bmax_ols srds_ols bmax srds bmax3 srds3],info)
disp(' ')
disp('The closest neighbor approach is intermediate to a non-spatial approach (OLS) and a full spatial approach (MESS)');


tmat=[closest_ar_time;MESS_time;ols_time];
disp(' ' )
disp('The following table shows the timings for the closed-form estimators -- all are very fast for 3,107 observations.')
disp(' ')

%making a table using James LeSage's mprint routine at
%www.spatial-econometrics.com
%info1.cnames=strvcat('+/- dependence','+ dependence');
info1.rnames=strvcat('ESTIMATION TIME TABLE: ', 'Closest AR', 'MESS', 'OLS');
mprint(tmat,info1)



