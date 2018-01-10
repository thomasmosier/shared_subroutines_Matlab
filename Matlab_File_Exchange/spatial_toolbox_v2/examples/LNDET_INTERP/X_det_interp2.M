%This example program examines timings from symmetric, asymmetric weight matrices using either positive dependence
%or positive and negative dependence.
%
%Written by Kelley Pace, www.spatial-statistics.com, 12/25/02

clear all;
clear global;
close all;
format short;
clc

%number of observations
n=1000;

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
[abound]=fmin_eig2(d);%generate minimum eigenvalue
ninterp=20;%number of interpolation points
ntotal=995;%ensures increments of 0.001 in a
[alphacoarse, alphafine]=fgrid_generate2(ninterp, ntotal, abound);%generates various grids


%Calling exact symmetric log-determinant routine for positive and negative spatial dependence
sym=1;
tic;
[detvalz]=fdet_interp2(d, sym, alphacoarse, alphafine);
s_posneg_time=toc;


%Calling exact asymmetric log-determinant routine for positive and negative spatial dependence
tic;
sym=0;
[detvalz]=fdet_interp2(d, sym, alphacoarse, alphafine);
a_posneg_time=toc;


%Calling exact symmetric log-determinant routine for positive spatial dependence
tic;
sym=1;
[detvalz]=fdet_interp2(d, sym);
s_pos_time=toc;


%Calling exact asymmetric log-determinant routine for positive spatial dependence
tic;
sym=0;
[detvalz]=fdet_interp2(d, sym);
a_pos_time=toc;

%putting times into a matrix
tmat(1,:)=[s_posneg_time s_pos_time];
tmat(2,:)=[a_posneg_time a_pos_time];

disp(' ' )
disp('The following table shows the timing penalties for adding negative dependence and assymetry')
disp(' ')
%making a table using James LeSage's mprint routine at
%www.spatial-econometrics.com
info1.cnames=strvcat('+/- dependence','+ dependence');
info1.rnames=strvcat('DETERMINANT TIME TABLE: ','symmetric','asymmetric');
mprint(tmat,info1)
disp(' ')
disp('For smaller n, the differences are not large in absolute terms. As n becomes large, it pays to') 
disp('take advantage of symmetry and the direction of spatial dependence (when known).')
