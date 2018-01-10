function [detvalz, alphafine]=fdet_interp_seq2(d, sym, alphacoarse, alphafine)
%
%[detvalz, alphafine]=fdet_interp_seq2(d, sym, alphacoarse, alphafine)
%
%This function computes ln|I-a*d(1:i,1:i)| for i=1 to n for a grid of values of a (stored in alphafine) where d is an
%n by n spatial weight matrix.
%
%
%INPUT:
%
%The n by n spatial weight matrix d is a required input.
%
%The optional scalar sym gives whether a matrix is asymmetric (sym=0) or symmetric (sym=1). Not supplying this forces
%the routine to test for symmetry which could cost time for large problems.
%
%The optional itersub by 1 vector alphacoarse gives the actual evaluation points of the log-determinant for the 
%corresponding values of a in alphacoarse. If you supply alphacoarse, you must supply alphafine.
%
%The optional aiter by 1 vector alphafine gives the arguments associated with the output log-determinant values. 
%If you supply alphafine, you must supply alphacoarse.
%
%
%OUTPUT:
%
%detvalz is a iter by n matrix where the first column is a vector of the
%values of the spatial dependence parameter (equal to alphafine when input), and the second column is a
%vector of the log-determinants.
%
%The aiter by 1 vector alphafine gives the arguments associated with the output log-determinant values. 
%
%
%NOTES:
%
%This function was discussed in:
%
%Pace, R. Kelley, and James LeSage, “Spatial Autoregressive Local Estimation,” 
%Spatial Statistics and Spatial Econometrics, Edited by Art Getis, Palgrave, 2003.
%
%
%Intially written by Kelley Pace, www.spatial-statistics.com, on 3/19/98, and revised 12/25/02.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% Default setting %%%%%%%%%%%%%%%%%%%%%%%

if (nargin<2)
    adev=max(max(abs(d-d')));
    if adev<1000*eps
        sym=1;
    else
        sym=0
    end
end


if (nargin<3)
%arbitrarily examine iter values of the log-determinant. This could be
%changed. This gives how many values of the log-determinant come from the
%function.

%since the output is iter by n, increasing iter to high values could result in large matrices
iter=100;
%This linearly spaces them between 0.01 and 0.99
alphafine=linspace(0.00, 0.99, iter)';
%The log-determinant is nonlinear, and so this spaces evaluation points equally for the function (appproximate).
alphacoarse=(1-exp(linspace(0, log(0.01), 21)))';
end

if (nargin==3)
    error('If you supply alphacoarse, you must supply alphafine')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%finding the size of the weight matrix
[n,ncols]=size(d);

%%%%%%%%%%%%%%%%% block of partial input checking %%%%%%%%%%%%%%%
%finding number of output points
aiter=length(alphafine);
%making sure alphafine is a column vector
      if aiter==1
       error('alphafine needs to be a column vector') 
    end
    
%transforming a row into a column vector (when necessary)   
      if size(alphafine,2)>1
        alphafine=alphafine';
    end
        
   if max(alphafine)>=1
       error('alphafine should have a maximum element of less than 1 and a minimum element of greater than the inverse of the minimum eigenvalue')
   end
   
%finding number of evaluation points
itersub=length(alphacoarse);
%making sure alphafine is a column vector
      if itersub==1
       error('alphacoarse needs to be a column vector') 
    end
    
%transforming a row into a column vector (when necessary)   
      if size(alphacoarse,2)>1
        alphacoarse=alphacoarse';
    end
        
   if max(alphacoarse)>=1
       error('alphacoarse should have a maximum element of less than 1 and a minimum element of greater than the inverse of the minimum eigenvalue')
   end
   
   if n~=ncols
     error('weight matrix needs to be n by n')
   end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%%%%%%%%%%%%%%%%% Preliminary Setup %%%%%%%%%%%%%%%%%%%%%%%

%these settings help optimize the sparse matrix routines
spparms('tight');
spparms('autommd',0);

%itersub gives how many interpolation points are used
itersub=length(alphacoarse);

%This forms a sparse identity matrix of order n
s1=speye(n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Asymmetric matrix %%%%%%%%%%%%%%%%%%%%%
if sym==0

%intitialization of interpolation point loop
detsub=zeros(itersub,n);
for i=1:itersub;
alpha=alphacoarse(i);

%LU decomposition of spatial transformation for asymmetric d.
[l,u]=lu(s1-alpha*d);
%log det is the sum of the log of the pivots
detsub(i,:)=cumsum(log(diag(u)))';

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Symmetric matrix %%%%%%%%%%%%%%%%%%%%%
if sym==1

%intitialization of interpolation point loop
detsub=zeros(itersub,n);
for i=1:itersub;
alpha=alphacoarse(i);

%Cholesky decomposition of spatial transformation for symmetric d.
u=chol(s1-alpha*d);
%In the determinant computations, the cholesky takes the pivots to the 1/2 power. Hence we multiply by 
%2 after taking the sum of the log of the pivots to yield the
%log-determinant.
detsub(i,:)=2*cumsum(log(diag(u)))';   

end

end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%interpolating for a finer grid of alpha
detvalz=spline(alphacoarse', detsub', alphafine')';

