function dm=fts(z, ms, rho, nnmat)
%
%This function computes the spatial lag of the variables z based on the index matrix nnmat.
%It gives a weighting directly related to the parameter rho.
%
%INPUT:
%
%The n by p matrix z represents the underlying set of variables
%
%The n by m matrix nnmat contains the indices of the neighboring observations
%
%The scalar ms represents the number of neighbors to use, and this should be greater than or equal to m.
%
%The scalar rho is a weighting parameter, where larger values lead to less decline in the weights with higher-order
%neighbors. It is 1-geometric decay rate.
%
%
%OUTPUT:
%
%The n by p matrix dm is the spatially averaged version of z. In other words, d*z where d is a spatial weight matrix.
%
%
%NOTES:
%
%This was used in:
%
%Pace, R. Kelley and Ronald Barry, O.W. Gilley, C.F. Sirmans, 
%“A Method for Spatial-temporal Forecasting with an Application to Real Estate Prices,” 
% International Journal of Forecasting, Volume 16, Number 2, April-June 2000, p. 229-246.
%
%Written by Kelley Pace, www.spatial-statistics.com, 1/9/03.

%computing weights given to neighbors
lagcoef=rho.^(1:ms);
%standardizing these to sum to 1
lagcoef=lagcoef/sum(lagcoef);

%size of input matrix
[n,k]=size(z);
%preallocating output matrix
dm=zeros(n,k);
%weighting and summing to form spatial average
for i=1:ms
    dm=z(nnmat(:,i), : )*lagcoef(i)+dm;
end;

    