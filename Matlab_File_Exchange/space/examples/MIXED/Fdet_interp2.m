function [detvalz]=fdet_interp2(d, sym, alphacoarse, alphafine)
%
% [detvalz]=fdet_interp2(d, sym, alphacoarse, alphafine)
%
%This function computes a vector of log-determinants for a vector of AR parameters (alpha)
%It uses a spline interpolation routine to reduce the number of determinants computed.
%
%
%INPUT:
%
%The n by n spatial weight matrix d is a required input.
%
%The optional scalar sym gives whether a matrix is asymmetric (sym=0) or symmetric (sym=1). Not supplying this forces
%the routine to test for symmetry which could cost time for large problems.
%
%The optional ncoarse by 1 vector alphacoarse gives the actual evaluation points of the log-determinant for the 
%corresponding values of a in alphacoarse. If you supply alphacoarse, you must supply alphafine.
%
%The optional nfine by 1 vector alphafine gives the output log-determinant values for the corresponding
%values of a in alphafine. If you supply alphafine, you must supply alphacoarse.
%
%
%OUTPUT:
%
%detvalz is a iter by 2 matrix where the first column is a vector of the
%values of the spatial dependence parameter (equal to alphafine when input), and the second column is a
%vector of the log-determinants.
%
%NOTES:
%
%The use of a grid of log-determinants, sparse weight matrices, and
%interpolation was discussed in:
%
%Pace, R. Kelley, and Ronald Barry, Quick Computation of Regressions with a Spatially Autoregressive Dependent Variable,
%Geographical Analysis, Volume 29, Number 3, July 1997, p. 232-247. 
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
iter=1000;
%This linearly spaces them between 0.01 and 0.99
alphafine=linspace(0.000, 0.999, iter)';
%The log-determinant is nonlinear, and so this spaces evaluation points equally for the function (appproximate).
alphacoarse=(1-exp(linspace(0, log(0.01), 21)))';
end

if (nargin==3)
    error('If you supply alphacoarse, you must supply alphafine')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
%We form a trial spatial transformation z=I-aD to aid in finding a
%ordering
z=s1-.1*d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Asymmetric matrix %%%%%%%%%%%%%%%%%%%%%
if sym==0
%We find an approximate column minimum degree ordering. One could use other
%permutations as indicated by the commented out commands below.
p=colamd(z);%approximate column minimum degree
%p=colmmd(z); %exact minimum degree

%intitialization of interpolation point loop
detsub=zeros(itersub,1);
for i=1:itersub;
alpha=alphacoarse(i);

%LU decomposition of spatial transformation.
z=s1-alpha*d;

%LU decomposition of spatial transformation for asymmetric d.
[l,u]=lu(z(:,p));
%log det is the sum of the log of the pivots
detsub(i)=sum(log(diag(u)));

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Symmetric matrix %%%%%%%%%%%%%%%%%%%%%
if sym==1

%We find an approximate minimum degree ordering. One could use other
%permutations as indicated by the commented out commands below.
ps=symamd(z);%approximate minimum degree
d=d(ps,ps);
%ps=symmmd(z); %exact minimum degree
%ps=symrcm(z); %reverse Cuthill-McKee

%intitialization of interpolation point loop
detsub=zeros(itersub,1);
for i=1:itersub;
alpha=alphacoarse(i);

%Cholesky decomposition of spatial transformation for symmetric d.
u=chol(s1-alpha*d);
%In the determinant computations, the cholesky takes the pivots to the 1/2 power. Hence we multiply by 
%2 after taking the sum of the log of the pivots to yield the
%log-determinant.
detsub(i)=2*sum(log(diag(u)));   
end

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%interpolating for a finer grid of alpha
detvalz=[alphafine spline(alphacoarse, detsub, alphafine) ];



%interpolating for a finer grid of alpha
detvalz=[alphafine spline(alphacoarse, detsub, alphafine) ];

