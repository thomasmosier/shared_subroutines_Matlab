%This example program uses the Chebyshev log-determinant approximation
%along with Taylor series bounds on the log-determinant. It also compares
%these to an exact log-determinant routine.
%
%Written by Kelley Pace, www.spatial-statistics.com, 12/25/02
%
%Takes 8 seconds to run on a 1 ghz PIII

clear all;
clear global;
close all;
format short;
clc

%number of observations
n=5000;

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

d=wswdel;
sym=1;
abound=-1;
ninterp=20;%number of interpolation points
ntotal=995;%ensures increments of 0.001 in a
[alphacoarse, alphafine]=fgrid_generate2(ninterp, ntotal, abound);%generates various grids

%Calling Chebyshev log-det approximation and Taylor series bounds
tic;
[detbounds]=fdet_chebyshev2(d,  alphafine);
chebyshev_time=toc;

%Calling exact routine
tic;
[detvalz]=fdet_interp2(d, sym, alphacoarse, alphafine);
exact_time=toc;

disp(' ')
disp('Chebyshev approximation time, exact log-determinant time')
[ chebyshev_time  exact_time]
disp(' ')
disp('Both the approximate and exact times are fast for this problem. However, for larger problems the');
disp('exact method will begin to rise almost qudratically with n while the Chebyshev will rise linearly with n.');
disp('Thus, for larger problems the Cheybshev will begin to show its strength.');
disp(' ');


%comparison in levels
plot(detbounds(:,1),detbounds(:,2),'--',detbounds(:,1),detbounds(:,3),'o',detbounds(:,1),detbounds(:,4),'+',detvalz(:,1),detvalz(:,2),'x')
xlabel('\alpha')
ylabel('log-determinants')
title('Exact log-determinant, Chebyshev approximation, Taylor bounds')
legend('Taylor Lower -', 'Chebyshev o', 'Taylor Upper +', 'Exact',0)


%comparison in differences
dcheby=diff(detbounds(:,3));
dexact=diff(detvalz(:,2));
ndseq=detbounds(2:end,1);

figure
plot(ndseq,dcheby,'o', ndseq,dexact,'x')
xlabel('\alpha')
ylabel('Difference in log-determinants')
title('Differences in exact log-determinant, Chebyshev approximation for successive \alpha')
legend('Chebyshev Differences',  'Exact Differences')

disp('The Chebyshev approximation appears quite good for positive values of the dependence parameter,');
disp('but could use improvement for materially negative values of the spatial dependence parameter.');
disp('Fortunately, such negative values seem rare in practice.');

