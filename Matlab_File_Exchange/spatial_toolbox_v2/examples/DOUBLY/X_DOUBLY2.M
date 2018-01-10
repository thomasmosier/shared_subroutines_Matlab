%This example program examines the doubly stochastic scaling routine.
%Written by Kelley Pace, www.spatial-statistics.com, 12/25/02

clear all;
clear global;
close all;
format short;
clc

%sample size
n=1000;

%unit vector
o=ones(n,1);

%creating random locational coordinates
xcoord=rand(n,1);
ycoord=rand(n,1);


%wswdel is a symmetric weight matrix
%wwsdel is a row-stochastic assymetric weight matrix formed from a symmetric adjacency matrix
%wmatdel is the diagonal weighting matrix used in making the above

[wswdel,wwsdel,wmatdel]=fdelw2(xcoord,ycoord);


tic
doubly_d=fdoubly2(wswdel);
douby_d_time=toc

disp('multiply douby_d times the unit vector o, which should approximately equal o. Here is the min, max doubly_d*o')
disp(' ')
[min(doubly_d*o) max(doubly_d*o)]

disp('this applies for other powers as well. Here is the second power.')
disp(' ')
[min(doubly_d*doubly_d*o) max(doubly_d*doubly_d*o)]

disp('min, max of an ordinary delaunay weight matrix times a unit vector')
disp(' ')
[min(wswdel*o) max(wswdel*o)]
disp(' ')
disp('Note, the material difference between the min and the max.');

