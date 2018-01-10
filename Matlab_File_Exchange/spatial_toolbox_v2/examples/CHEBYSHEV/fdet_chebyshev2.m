function [detbounds]=fdet_chebyshev2(d, alphafine)
%
%[detbounds]=fdet_chebyshev2(d, alphafine)
%
%This function returns a quadratic lower Taylor series bound, a quadratic Chebyshev approximation, and an upper
%Taylor series bound for the log-determinant of (I-aD).  As a default the elements of a can go from 0 to 1 or can take the optional
%vector alphafine containing elements between 1 and the reciprocal of the minimum eigenvalue of d. The computation of 
%log-determinants for negative values of a does not require more time, but may cause other routines to take more time.
%
%
%INPUT:
%
%The n by n symmetric weight matrix d with zeros on the diagonal.
%
%The optional n by 1 vector alphafine containing values for a. These values should lie between 
%1 and the reciprocal of the minimum eigenvalue of d. As a default, alphafine has 1,000 elements with
%the initial element at 0 and the largest element is 0.999.
%
%
%OUTPUT:
%
%detbounds is a iter by 4 matrix comprised of:
%a iter by 1 grid of values of the dependence parameter a.
%a iter by 1 vector of the lower Taylor bound for the log-determinant
%a iter by 1 vector of the quadratic Chebyshev approximation
%a iter by 1 vector of the upper Taylor bound for the log-determinant
%
%NOTES:
%
%If you use this, please cite:
%
%Pace, R. Kelley, and James LeSage, 
%"Chebyshev Approximation of Log-determinants of Spatial Weight Matrices,"
%Computational Statistics and Data Analysis, forthcoming.
%
%and
%
%Pace, R. Kelley, and James P. LeSage, 
%"Semiparametric Maximum Likelihood Estimates of Spatial Dependence," 
%Geographical Analysis, Volume 34, Number 1, January 2002, p. 75-90.
%
%Kelley Pace, www.spatial-statistics.com, December 7, 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Default setting %%%%%%%%%%%%%%%%%%%%

if (nargin<2)
%alpha vector default
alphafine=(0.0:0.001:0.999)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%finding the size of the weight matrix
[n,ncols]=size(d);
aiter=length(alphafine);

%%%%%%%%%%%%%%%%% block of partial input checking %%%%%%%%%%%%%%%

%making sure alphafine is a column vector
      if aiter==1
       error('alphafine needs to be a column vector') 
    end
    
%transforming a row into a column vector (when necessary)   
      if size(alphafine,2)>1
        alphafine=alphafine';
    end
    
    if n~=ncols
     error('weight matrix needs to be n by n')
   end
   
   if max(alphafine)>=1
       error('alphafine should have a maximum element of less than 1 and a minimum element of greater than the inverse of the minimum eigenvalue')
   end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%td1=trace(d); If one needs to compute this. Usually is zero by
%construction. The Taylor bounds assume this.
td1=0;

%trace(D*D)
%For really large problems, one might take only one of the triangles, square elements, sum, and scale by 2
td2=sum(sum(d.^2));

%first row contains the coefficients of T(0)
%second row contains the coefficients of T(1)
%third row contains the coefficients of T(2)
%where T(0), T(1), T(2) represent the first three Chebyshev polynomials
rezzy=[[1 0 0];[0 1 0];[-1 0 2]];

%This is 1 plus the highest order (hence, quadratic=3)
nposs=3;
%A column vector sequence 1,2, and 3
seq1nposs=[1;2;3];

%The following equations compute a discrete Chebshev approximation for each value of the spatial
%dependence vector rho. This follows the development in Press et al.,
%Numerical Recipes in Fortran 77, second edition, Cambridge University
%Press and applies it to the function sum(ln(1-rho*xk)) 
xk=cos(pi*(seq1nposs-0.5)/nposs);

for i=1:aiter; 
rho=alphafine(i);   
for j=1:nposs;
cposs(i,j)=(2/nposs)*sum(log(1-rho*xk).*cos(pi*(j-1).*(seq1nposs-0.5)/nposs));
end;
end;

%This multiplies the Chebyshev coefficients times the polynomials
comboterm=cposs*rezzy;

%This is the size, trace(D), trace(D*D)-n/2] -- for most spatial weight
%matrices trace(D)=td1=0. The Taylor bounds assume this. The last element
%is also adjusted by -0.5*n which appears in the approximation formula.
tdvec=[n;td1;(td2-0.5*n)];

%This computes the actual Chebyshev quadratic approximation of the
%log-determinant for each value of rho in the grid
lnest=comboterm*tdvec;

%Taylor series bounds
%altlowquad=(log(1-alphafine)+alphafine)*td2;
%altupquad=-(alphafine.^2)*0.5*td2;

posind=(alphafine>=0);

altlowquadpos=(log(1-alphafine)+alphafine)*td2;
altupquadpos=-(alphafine.^2)*0.5*td2;

altlowquadneg=log(1-alphafine.^2)*0.5*td2;
altupquadneg=(log(1-alphafine)+alphafine)*td2;

altlowquad=posind.*altlowquadpos+(1-posind).*altlowquadneg;
altupquad=posind.*altupquadpos+(1-posind).*altupquadneg;

%Both Chebyshev Approximation and Taylor bounds
detbounds=full([alphafine altlowquad lnest altupquad]);



