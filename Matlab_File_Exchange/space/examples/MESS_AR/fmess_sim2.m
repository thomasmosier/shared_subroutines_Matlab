function [rvcorr]=fmess_sim2(d,rv,alpha)
%
%function [rvcorr]=fmess_sim2(d,rv,alpha)
%
%This function simulates a normal random variable with a matrix exponential variance-covariance function
%In other words, it simulates N(0,expm(alpha*d)), where expm is the matrix exponential function.
%
%INPUT:
%
%d is a n by n spatial weight matrix. The function assumes this is
%symmetric.
%
%rv is a n by iter matrix of random deviates
%
%alpha is the spatial dependence parameter
%
%OUTPUT:
%
%rvcorr is the n by iter matrix of spatially correlated random variates
%
%NOTES:
%
%If you use the function, please cite:
%
%LeSage, James and R. Kelley Pace, Spatial Dependence in Data Mining,
%Data Mining for Scientific and Engineering Applications, 
%Edited by Robert L. Grossman, Chandrika Kamath, Philip Kegelmeyer, Vipin Kumar, and Raju R. Namburu, 
%Kluwer Academic Publishing, 2001.
%
%Kelley Pace, www.spatial-statistics.com, 12/25/02

%In the matrix exponential world one can take the square root of the matrix exponential by simply
%taking the matrix exponential with one-half the parameter of the original matrix.
truerho=alpha/2;

%maximum order of D to use in computing the matrix exponential
%a sufficient condition of weight matrices with maximum eigenvalue of 1 is
%to make (max(abs(alphavec))^(q))/(q)! a small number (eg. 0.005)
%lower this to reduce time, if the max(abs(alphavec)) is sufficiently low
q=24;

%sequence of powers of the spatial dependence parameter from 0 to q-1
nr=(truerho.^(0:(q-1)))';
%sequence of associated factorials from 0:q-1
dr=gamma(1:q)';
%MESS weights
rr=nr./dr;

%initialization where i=1
rvcum=rv;
rvcorr=rv;

for i=2:q
%gives (D^i)*rv    
rvcum=d*rvcum;	
%sums weighted powers of D times the input random variates	
rvcorr=rvcum*rr(i)+rvcorr;
end;



	