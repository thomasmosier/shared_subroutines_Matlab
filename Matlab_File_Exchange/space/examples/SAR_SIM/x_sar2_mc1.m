%This simulates and estimates a simultaneous autoregression (sar)
%
%Written by Kelley Pace, www.spatial-statistics.com, 12/24/02

clear all;
clear global;
clc;
close all;

n=10000;
p=1;
k=p+1;
iter=1;
alpha=0.8;
sig=1;
info.cnames=strvcat('beta est','SRDS   ','  Probs of SRD+');
info.rnames=strvcat(' ','Beta 1','Intercept','Alpha');

%Truth
disp(' ')
truevec=[n; k ;1 ;alpha];
infot.rnames=strvcat('Parameters, Variables','n','k','beta','alpha');
mprint(truevec,infot);

%Simlated variables
xcoord=rand(n,1);
ycoord=rand(n,1);
rv=randn(n,iter);
xsub=randn(n,p);
o=ones(n,1);
x=[xsub o];
beta=ones(k,1);

%finding Delaunay weight matrices
tic;
[wswdel,wwsdel,wmatdel]=fdelw2(xcoord,ycoord);
weight_matrix_time=toc;

%selecting a symmetric weight matrix
d=wswdel;

%Correlated deviates
tic;
[rvcorr]=fsar_sim2(d,rv,alpha);
sim_time=toc;

%simulated dependent variable
y=x*beta+sig*rvcorr;



%Estimating a sar model without log-determinant 
%using Chebyshev approximations, Taylor bounds, and likelihood dominance inference
tic
 [bmax1, srds1, prhighers1, emax1, logliks1]=fsar2(x, y, d);
time_sans_lndet=toc;

ests1=[bmax1 srds1 prhighers1];

disp('fsar2 Estimation Results using approximate lndet');
mprint(ests1,info)
disp(' ')


%computing log-determinants
tic;
[detvalz]=fdet_interp2(d);
log_det_time=toc;
 

%Estimating sar
tic
 [bmax, srds, prhighers, emax, logliks]=fsar2(x, y, d, detvalz);
time_with_lndet=toc;

ests=[bmax srds prhighers];

disp('fsar2 Estimation Results using exact lndet');
mprint(ests,info)
disp(' ')
disp('The approximate approach shows similar estimates to the exact approach.');
disp('Note, the likelihood dominance SRDs are smaller in magnitude than the exact SRDs.');
disp('However, they can still document statistical significance for many variables, and thus can');
disp('prove useful in many circumstances.');
disp(' ')
disp(' ')


xmix=[xsub d*xsub o];
tic;
[bmax2, srds2, prhighers2, emax2, logliks2]=fmix2(xmix, y, d, detvalz);
mix_time=toc;

estmix=[bmax2 srds2 prhighers2];
infomix.cnames=strvcat('beta est','SRDS   ','  Probs of SRD+');
infomix.rnames=strvcat(' ','Beta for X','Beta for DX','Intercept','Alpha');

disp('The mixed model subsumes SAR and should show similar alpha estimates.');
disp('In addtion, the regression estimate for DX should equal -alpha if the DGM is SAR.');
mprint(estmix,infomix)
disp(' ')



disp('Times (must add lndet time to exact sar to compare to approximate sar)');
time_vec=[weight_matrix_time;sim_time;log_det_time;time_sans_lndet;time_with_lndet;mix_time];
info1.rnames=strvcat('Times','Weight Matrix','Simulation','Lndet','Approximate sar','Exact sar','mixed time');
mprint(time_vec,info1)






 






 

 
   
 
 
 




