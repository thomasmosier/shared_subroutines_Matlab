function [wwsasymnn,wmatasymnn]=fasym_neighbors2(m,rho)
%
%[wwsasymnn,wmatasymnn]=fasym_neighbors2(m,rho)
%
%This function takes the individual neighbor matrices and manipulates
%these to form an assymmetric, normalized weight matrix based on m nearest
%neighbors weighted according to rho.
%
%INPUT:
%
%The scalar m which gives how many neighbors to use
%
%The scalar rho which sets the geometric weighting (1-rho)=decay
%
%The file smats.mat on disk from a previous invocation of fneighbors2.
%
%OUTPUT:
%
%wswnn is a symmetric n by n spatial matrix similar to (same eigenvalues) the asymmetric row-stochastic wwsnn
%
%wwsasymnn is a asymmetric row-stochastic weight matrix
%
%wmatasymnn is the n by n diagonal matrix so that wmatasymnn*a is row-stochastic.
%
%NOTES:
%
%One can optimize the likelihood over neighbors and weighting of neighbors
%as in:
%
%Pace, R. Kelley and Ronald Barry, O.W. Gilley, C.F. Sirmans, 
%“A Method for Spatial-temporal Forecasting with an Application to Real Estate Prices,” 
%International Journal of Forecasting, Volume 16, Number 2, April-June 2000, p. 229-246.
%
%Written by Kelley Pace, www.spatial-statistics.com, on 6/23/97 and revised
%12/26/02


%Loading from disk the individual neighbor matrices created by the function fneighbors2
load smats;

[n,n]=size(s1);
sall=sparse(n,n);
%rho is 1 minus the rate of decay of weights with distance
%m is the number of neighbors
%this creates the weights
w=rho.^(0:(m-1));
%this block accumulates the weighted, asymmetric matrices
for i=1:m;
   si=eval(['s' num2str(i)]);
   sit=(si>0);
   sall=w(i)*sit+sall;
end;
%computes the row sum
sumvec=sum(sall')';
%creates normalizing diagonal matrix
wmatasymnn=spdiags((1./sumvec),0,n,n);
%creates row-stochastic matrix with 1=max(eigenvalue)
wwsasymnn=wmatasymnn*sall;
