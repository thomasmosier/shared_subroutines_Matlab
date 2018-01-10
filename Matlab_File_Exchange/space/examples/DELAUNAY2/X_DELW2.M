%This example program uses the Delaunay triangularization routine in Matlab to arrive at a spatial weight matrix.
%
%This example program uses the data from:
%Pace, R. Kelley, and Ronald Barry, Quick Computation of Regressions with a Spatially Autoregressive Dependent Variable,
%Geographical Analysis, Volume 29, Number 3, July 1997, p. 232-247. 
%
%Written by Kelley Pace, www.spatial-statistics.com, 6/23/97 revised 12/25/02

clear all;
clear global;
close all;
format short;
clc

%locations on the plane
load xcoord;
load ycoord;

%number of observations
n=length(ycoord);

%wswdel is a symmetric weight matrix
%wwsdel is a row-stochastic assymetric weight matrix formed from a symmetric adjacency matrix
%wmatdel is the diagonal weighting matrix used in making the above
tic;
[wswdel,wwsdel,wmatdel]=fdelw2(xcoord,ycoord);
delaunay_time=toc

%This plots the graph (as in graph theory -- edges and nodes) underlying the weight matrix 
gplot(wswdel,[xcoord ycoord])
xlabel('Longitude in Decimal Degrees')
ylabel('Latitutde in Decimal Degrees')
title('Delaunay Nearest Neighbor Graph')

%One can visualize the weight matrix -- the patches arise due to the
%grouping of counties by state. 
figure
spy(wswdel)
xlabel('Observation Number')
ylabel('Observation Number')
title('Original Non-zero Weight Pattern for Delaunay')

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
title('Permuted Non-zero Weight Pattern for Delaunay')


