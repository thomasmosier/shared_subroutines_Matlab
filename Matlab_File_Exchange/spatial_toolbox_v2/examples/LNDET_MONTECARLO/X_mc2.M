%This example program uses the Monte Carlo log-determinant approximation
%with confidence limits. It also compares
%these to an exact log-determinant routine.
%
%Written by Kelley Pace, www.spatial-statistics.com, 1/03/03

clear all;
clear global;
close all;
format bank;
clc

%number of observations
n=2500;

%random points on a plane
xcoord=rand(n,1);
ycoord=rand(n,1);

%Creating spatial weight matrices
%wswdel is a symmetric weight matrix
%wwsdel is a row-stochastic assymetric weight matrix formed from a symmetric adjacency matrix
%wmatdel is the diagonal weighting matrix used in making the above
tic;
[wswdel,wwsdel,wmatdel]=fdelw2(xcoord,ycoord);
delaunay_time=toc;

%[mineig]=fmin_eig2(wswdel);%generate minimum eigenvalue
%abound=1/mineig+0.005;

abound=-0.995;
ninterp=20;%number of interpolation points
ntotal=995;%insures increments of 0.001 in spatial dependence parameter
%generating both negative and positive spatial dependence for interpolation grid and output grid.
 [alphacoarse, alphafine]=fgrid_generate2(ninterp, ntotal, abound);

%Calling Monte Carlo log-det approximation with confidence limits.
tic;
%confide=fdet_mc2(wswdel);%default with positive dependence
confide=fdet_mc2(wswdel, 30, 16, 0.01, 1111, alphafine);%fully specified
mc_time=toc;

%Calling exact routine
tic;
%[detvalz]=fdet_interp2(wswdel);%default with positive dependence
[detvalz]=fdet_interp2(wswdel, 1, alphacoarse, alphafine);%fully specified
exact_time=toc;

disp(' ')
disp('Monte Carlo approximation time, exact log-determinant time')
[ mc_time  exact_time]
disp(' ')
disp('Both the approximate and exact times are fast for this problem. However, for larger problems the');
disp('exact method will begin to rise almost qudratically with n while the Monte Carlo will rise linearly with n.');
disp('Thus, for larger problems the Monte Carlo will begin to show its strength.');
disp(' ');


%comparison in levels
plot(confide(:,1),confide(:,2),'--',confide(:,1),confide(:,3),'o',confide(:,1),confide(:,4),'+',detvalz(:,1),detvalz(:,2),'x')
xlabel('\alpha')
ylabel('log-determinants')
title('Exact log-determinant, Monte Carlo approximation and confidence limits')
legend('Lower Confidence', 'Monte Carlo', 'Upper Confidence', 'Exact',0)


disp('The Monte Carlo approximation performs excellently!');
