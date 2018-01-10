function tnnmat=fnn_mex2(xcoord, ycoord, ms, n, istart)
%
%tnnmat=fnn_mex2(xcoord, ycoord, ms, n, istart)
%
%This function computes m spatiotemporal nearest neighbors for each observation from observation 1 to (n-istart). 
%This matlab function calls a Fortran 90 Mex function, tnmex4.dll. For other platforms, the fortran source code is
%available in a subdirectory within the space_time subdirectory under the examples subdirectory.
%
%INPUT:
%
%xcoord, ycoord are n by 1 vectors ordered by time with observation 1 the most recent and observation n the oldest.
%
%ms is the number of nearest neighbors previous in time
%
%n is the number of observations
%
%istart is where to begin computing temporal nearest neighbors. This allows istart observations to exist before beginning, 
%and istart>ms.
%
%
%OUTPUT:
%
%tnnmat is a n by ms matrix of spatiotemporal nearest neighbors.
%
%
%NOTES:
%
%
%A similar function was used in:
%
%Pace, R. Kelley and Ronald Barry, O.W. Gilley, C.F. Sirmans, 
%“A Method for Spatial-temporal Forecasting with an Application to Real Estate Prices,” 
% International Journal of Forecasting, Volume 16, Number 2, April-June 2000, p. 229-246.
%
%Written by Kelley Pace, www.spatial-statistics.com, 1/08/03.

%increase m by 1 to allow for non-own neighbors
ms=ms+1;
%call Fortran 90 MEX function
nnmat=tnmex4(xcoord,ycoord,ms,n,istart);
%transpose and eliminate index to own observation
nnmat=nnmat';
tnnmat=nnmat(:,2:end);
