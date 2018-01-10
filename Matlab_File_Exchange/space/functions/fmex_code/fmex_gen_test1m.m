clear all
clear global
close all
format bank
clc

%MEXing via Lahey LF95 Fortran 95 compiler

!lf95 tnmex4.f90 fmexdyn4.f90 @mexlink.rsp -out tnn13.dll


%example test variables
n=5000;
xcoord=rand(n,1);
ycoord=rand(n,1);
ms=29;
istart=300;

%testing
tic;
nnmat=tnmex4(xcoord,ycoord,ms,n,istart);
toc


