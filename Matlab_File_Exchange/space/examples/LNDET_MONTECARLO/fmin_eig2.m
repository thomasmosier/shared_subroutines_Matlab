function [mineig]=fmin_eig2(d)
%
%[mineig]=fmin_eig2(d)
%
%This function computes an estimate of the minimum eigenvalue (accuracy specified to 0.005). 
%
%INPUT:
%
%A symmetric n by n spatial weight matrix d.
%
%OUTPUT:
%
%The scalar minimum eigenvalue.
%
%NOTES:
%
%Written by Kelley Pace, www.spatial-statistics.com,12/31/02.

OPTS.tol=0.005;
 OPTS.disp=0;
 mineig=eigs(d,1,'SA',OPTS);