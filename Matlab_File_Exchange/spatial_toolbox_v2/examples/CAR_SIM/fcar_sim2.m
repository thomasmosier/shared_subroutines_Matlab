function [rvcorr]=fcar_sim2(d, rv, alpha)
%
%[rvcorr]=fcar_sim2(d, rv, truerho)
%
%This function creates correlated random variables, rvcorr, of the same dimensions as the
%supplied random variables, rv. The rvcorr are conditional autoregressive random variables.
%
%
%INPUTS:
%
%d is a n by n symmetric weight matrix
%
%rv is a n by iter matrix of random variables
%
%alpha is the true value of autoregressive parameter.
%
%OUTPUT:
%
%rvcorr is a n by iter matrix of correlated random variates
%
%
%NOTES:
%
%this function uses n*iter elements and thus has problems for large n and/or large iter
%You may wish to unroll this function and place it in a program if you are creating
%large numbers of correlated random variables. The approximate minimum degree algorithm may or may
%not be optimal for your problem (symamd). You can change this easily. Repeated use of 
%this function is inefficient since it creates the permutation and the cholesky 
%decomposition each time it is invoked. Storing these would be much more efficient for
%larger problems.
%
%Written by Kelley Pace, www.spatial-statistics.com, on 6/23/97, revised 1/1/03.


[n,iter]=size(rv);
s1=speye(n);
%creates n by n sparse identity matrix

z=s1-alpha*d;
%creates I-aD

p=symamd(z);
%finds symmetric minimum degree permutation

[ps,pinvind]=sort(p);
pinvind=(pinvind)';
%finds inverse permutation

r=chol(z(p,p));
%performs Cholesky decomposition on permuted I-aD

rvcorr=r\rv;
%solves for correlated random variables

rvcorr=rvcorr(pinvind,:);
%inverts permutation 


