function tnnmat=fnn_ml2(xcoord, ycoord, m, n, istart)
%
%fnn_ml2(xcoord, ycoord, m, n, istart)
%
%This function computes m spatiotemporal nearest neighbors for each observation from observation 1 to (n-istart).
%
%INPUT:
%
%xcoord, ycoord are n by 1 vectors ordered by time with observation 1 the most recent and observation n the oldest.
%
%m is the number of nearest neighbors previous in time
%
%n is the number of observations
%
%istart is where to begin computing temporal nearest neighbors. This allows istart observations to exist before beginning, 
%and istart>m.
%
%
%OUTPUT:
%
%tnnmat is a n by m matrix of spatiotemporal nearest neighbors.
%
%
%NOTES:
%
%There is a fortran mex version of a similar function, and it runs faster. This function will be acceptably fast for moderate sample sizes
%in version 6.5 due to the JIT compiler.
%
%This function requires more or n-squared memory locations, and hence will not work for large n.
%
%A similar function was used in:
%
%Pace, R. Kelley and Ronald Barry, O.W. Gilley, C.F. Sirmans, 
%“A Method for Spatial-temporal Forecasting with an Application to Real Estate Prices,” 
% International Journal of Forecasting, Volume 16, Number 2, April-June 2000, p. 229-246.
%
%Written by Kelley Pace, www.spatial-statistics.com, 1/08/03.

%preallocate
d2mat=zeros(n,n);

%Ironically, faster due to pure scalar functions in 6.5. Computes squared-distances to all previous observations
for j=1:n
    for i=(j+1):n
        d2mat(i,j)=(xcoord(i)-xcoord(j))^2+(ycoord(i)-ycoord(j))^2;
    end
end

%preallocates
tnnmat=ones(n,m);

%finds nearest neighbors
for i=istart:n
[d2vecs,nnind]=sort(d2mat(i,1:(i-1)));
tnnmat(i,:)=nnind(1:m);
end
