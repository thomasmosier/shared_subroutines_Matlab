function [lowerbounds, chebyshevest, upperbounds, alphafine]=fdet_chebyshev_seq2(d, alphafine)
%
%[lowerbounds, chebyshevest, upperbounds, alphafine]=fdet_chebyshev_seq2(d, alphafine)
%
%%This function returns a sequence of quadratic lower Taylor series bounds, quadratic Chebyshev approximations, and upper
%Taylor series bounds for the log-determinant of (I-aD(1:i,1:i)) for i=1...n. The routine
%assumes the input matrix d is symmetric with a zeros on the diagonal.  As a default the elements of  a  go from 0 to 1 or can take the optional
%vector alphafine containing elements between 1 and -1. The computation of 
%log-determinants for negative values of a does not require more time, but may cause other routines to take more time.
%
%
%INPUT:
%
%The n by n symmetric weight matrix d with zeros on the diagonal.
%
%The optional n by 1 vector alphafine containing values for a. These values should lie between 
%1 and -1 (for this routine). As a default, alphafine has 1,00 elements with
%the initial element at 0 and the largest element is 0.99.
%
%
%OUTPUT:
%
%alphafine is a iter by 1 grid of values of the dependence parameter a ordered from negative to positive
%lowerbounds is a iter by n matrix of the lower Taylor bound for the log-determinant.
%chebyshevest is a iter by n matrix of the quadratic Chebyshev approximation.
%upperbounds is a iter by n matrix of the upper Taylor bound for the log-determinant.
%
%
%NOTES:
%
%The chebyshev approximation appeared in:
%
%Pace, R. Kelley, and James LeSage, 
%"Chebyshev Approximation of Log-determinants of Spatial Weight Matrices,"
%Computational Statistics and Data Analysis, forthcoming.
%
%The use of determinant sequences appeared in:
%
%Pace, R. Kelley, and James LeSage, “Spatial Autoregressive Local Estimation,” 
%Spatial Statistics and Spatial Econometrics, Edited by Art Getis, Palgrave, 2003.
%
%Kelley Pace, www.spatial-statistics.com, December 7, 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Default setting %%%%%%%%%%%%%%%%%%%%

if (nargin<2)
%alpha vector default -- shouldn't be too large since output is n times length of this vector
%if there are negative values, they should be first and followed by positive values
alphafine=(0.0:0.01:0.99)';
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
%we need a row unit vector
o=ones(1,n);

%we need a row sequence vector
nnseq=(1:n);

%td1=trace(d); If one needs to compute this. Usually is zero by
%construction. The Taylor bounds assume this.
td1=0;

%trace(D(1:i)*D(1:i)) for sequence i=1 ... n
d_upper_SSQ=2.0*sum(triu(d).^2);
td2_sequence=cumsum(d_upper_SSQ);

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

%This is a sequence of size, trace(D), trace(D*D)-n/2] 
%The last vector is also adjusted by -0.5*nnseq which appears in the approximation formula.
lastvec=td2_sequence-0.5*nnseq;
tdvec=[nnseq;td1*o;lastvec];

%This computes the actual Chebyshev quadratic approximation of the
%log-determinant for each value of rho in the grid and for the sequence
%1...nn. lnest should have size iter by nn
chebyshevest=comboterm*tdvec;

%Taylor series bounds
posind=(alphafine>=0);
negind=logical(1-posind);

altlowquadpos=(log(1-alphafine)+alphafine)*td2_sequence;
altupquadpos=-(alphafine.^2)*0.5*td2_sequence;

altlowquadneg=log(1-alphafine.^2)*0.5*td2_sequence;
altupquadneg=(log(1-alphafine)+alphafine)*td2_sequence;

lowerbounds=[altlowquadneg(negind,:);altlowquadpos(posind,:)];
upperbounds=[altupquadneg(negind,:); altupquadpos(posind,:)];






