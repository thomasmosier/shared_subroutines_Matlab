%This example program uses the data from:
%Pace, R. Kelley, and Ronald Barry, Quick Computation of Regressions with a Spatially Autoregressive Dependent Variable,
%Geographical Analysis, Volume 29, Number 3, July 1997, p. 232-247. 

%This program computes nearest neighbor spatial weight matrices based upon the 1st and 2nd order Delaunay triangle neighbors
%If one requests more neighbors than found in the first and second order neighbors for a particular observation, the neighbors beyond these will have
%zero elements in the individual neighbor weight matrices. Hence, this is a set of nearest neighbors restricted to the 1st and 2nd order Delaunay neighbors.

%Important parameters include the number of neighbors, m, and rho, the decline of the weight of the jth neighbor with j

%Written by Kelley Pace, www.spatial-statistics.com, 6/23/97 and revised
%12/26/02

clear all;
clear global;
close all;
format short;
clc

%loading election data
load y;
load xsub;
load xcoord;
load ycoord;

%number of observations
n=length(ycoord);

%Computes Delaunay triangle based weight matrix and records timings
tic;
[wswdel,wwsdel,wmatdel]=fdelw2(xcoord,ycoord);
delaunay_time=toc;

%m is the number of neighbors desired
%rho specifies the geometrically declining weights
nnmax=30;
rho=0.85;


%needs delaunay weights
%Computes nearest neighbors and records timings
delorder=4;
tic;
 fneighbors2(wswdel, xcoord, ycoord, nnmax, delorder)
neighbor_time=toc;


m=8;
%Symmetricizes and weights NNs and records timings
[wswnn,wwsnn,wmatnn]=fsym_neighbors2(m, rho);
%wswnn is a symmetric weight matrix
%wwsnn is a row-stochastic assymetric weight matrix
%wmatnn is the diagonal weighting matrix used in making the above

%This plots the graph (as in graph theory -- edges and nodes) underlying the weight matrix 
gplot(wswnn, [xcoord ycoord])
xlabel('Longitude in Decimal Degrees')
ylabel('Latitutde in Decimal Degrees')
title('8 Nearest Neighbor Graph')


%One can visualize the weight matrix -- the patches arise due to the
%grouping of counties by state. 
figure
spy(wswnn)
xlabel('Observation Number')
ylabel('Observation Number')
title('Original Non-zero Weight Pattern for 8 Nearest Neighbors')


%I add a diagonal to the weight matrix to mimic the usual spatial transformations.
z=speye(n)+wswdel;

%One can permute the rows and columns of the weight matrix for lower
%bandwidth. This does not affect the determinant, but notice the sparsity
%of the matrix.
p=symrcm(z);%symmetric reverse Cuthill-McKee permutation

figure
spy(wswdel(p,p))
xlabel('Observation Number')
ylabel('Observation Number')
title('Permuted Non-zero Weight Pattern for 8 Nearest Neighbors')


%we can see the effect of more neighbors
m=25;
%Symmetricizes and weights NNs and records timings
tic;
[wswnn,wwsnn,wmatnn]=fsym_neighbors2(m, rho);
sym_neighbor_time=toc;
%wswnn is a symmetric weight matrix
%wwsnn is a row-stochastic assymetric weight matrix
%wmatnn is the diagonal weighting matrix used in making the above

%This plots the graph (as in graph theory -- edges and nodes) underlying the weight matrix 
figure
gplot(wswnn, [xcoord ycoord])
xlabel('Longitude in Decimal Degrees')
ylabel('Latitutde in Decimal Degrees')
title('25 Nearest Neighbor Graph')


%One can visualize the weight matrix -- the patches arise due to the
%grouping of counties by state. 
figure
spy(wswnn)
xlabel('Observation Number')
ylabel('Observation Number')
title('Original Non-zero Weight Pattern for 25 Nearest Neighbors')

%I add a diagonal to the weight matrix to mimic the usual spatial transformations.
z=speye(n)+wswdel;

%One can permute the rows and columns of the weight matrix for lower
%bandwidth. This does not affect the determinant, but notice the sparsity
%of the matrix.
p=symrcm(z);%symmetric reverse Cuthill-McKee permutation

figure
spy(wswdel(p,p))
xlabel('Observation Number')
ylabel('Observation Number')
title('Permuted Non-zero Weight Pattern for 25 Nearest Neighbors')


time_mat=[delaunay_time;neighbor_time;sym_neighbor_time];

disp(' ')
info1.rnames=strvcat('Timings in seconds for 3,107 observations','delaunay','nearest neighbors','symmetricizing and scaling');
mprint(time_mat, info1);



