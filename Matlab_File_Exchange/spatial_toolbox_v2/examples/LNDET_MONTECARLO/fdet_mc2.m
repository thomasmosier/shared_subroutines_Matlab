function confide=fdet_mc2(d, total_order, iter, asize, randstate, alphafine)
%
%confide=fdet_mc2(d, total_order, iter, asize, randstate, alphafine)
%
%This function estimates the Log-determinant of (I-a*d) and provides confidence bounds on the estimated log-determinant
%(a is a scalar, and d is an n by n weight matrix). While the estimated log-determinant is valid for a in (-1/min(eig(d)),1), the 
%lower bound computation requires a in (-1,1).
%
%
%INPUT:
%
%The required n by n matrix d can be a symmetric or asymmetric weight matrix. However, maximum modulus of the eigenvalue of d must be 1 or less. 
%A row-stochastic, doubly stochastic scaling will do this. Matrices similar to stochastic matrices will satisfy this as well.
%Finally, one can find the maximum eigenvalue, and scale a candidate weight matrix by this to yield a weight
%maxtrix with maximum eigenvalue of 1. Note, one can find the maximum eigenvalue of a sparse d using the command 
%eigs in a reasonable time for moderately large matrices. 
%
%The optional scalar integer total_order specifies the highest term estimated in the Taylor series (e.g., E(tr(D^total_order)) ). 
%Please specify even order. The default equals 30.
%
%The optional scalar integer iter specifies how many independent realizations of u to use in the estimation of E(tr(D^total_order)).
%The default equals 16.
%
%The optional scalar asize gives the size of the confidence limit desired in decimal units. Default is 0.05 (5%). This variable asize
%lie in (0,1). This is a Chebyshev interval, and the actual size is very conservative.
%
%The optional non-negative integer scalar randstate sets the state of the random number generator. Use this if you wish the same
%result each invocation, otherwise the results will vary somewhat each time. This variation should provide some idea of the 
%sensitivity of the results to the approximation.
%
%The optional aiter by 1 vector alphafine gives the evaluation points a(i) for i=1...aiter. The default is an linearly spaced vector
%going from 0 to 0.999 (aiter=1,000). It costs nothing to obtain log-determinant estimates for negative values of a in this routine,
%but negative values of a may lead other routines to perform more work. At the moment, the lower bound computation only works for
%a in the interval (-1,1).
%
%
%OUTPUT:
%
%The aiter by 4 matrix confide gives a, lower bound of ln|I-a*d|, an estimated ln|I-a*d|, and upper bound of ln|I-a*d|, 
%where aiter is specified internally (can be edited in this file) or has the same length of the
%optionally specified alphafine. Essentially, the matrix confide supplies the points used in the creation of the grid of 
%stimated log-determinant values along with their confidence limits. 
%
%confide(:,1) gives the evaluation points (values of a).
%
%confide(:,2) gives the lower confidence limit on the log-determinant
%
%confide(:,2) gives the estimated log-determinant
%
%confide(:,4) gives the upper confidence limit on the log-determinant
%
%
%NOTES:
%
%If you use this function, please cite:
%
%Barry, Ronald, and R. Kelley Pace, 
%"A Monte Carlo Estimator of the Log Determinant of Large Sparse Matrices," 
%Linear Algebra and its Applications, Volume 289, Number 1-3, 1999, p. 41-54.
%
%Written by Kelley Pace, www.spatial-statistics.com, on 11/27/99, revised 1/1/03.


%%%%%%%%%%%%%% Default setting %%%%%%%%%%%%%

if (nargin<6)
%alpha vector default
alphafine=(0.0:0.001:0.999)';
end

if (nargin>4)
rand('state', randstate);
end

if nargin<4
    %default size of test
    asize=0.05;
end

if nargin<3
%iter default
iter=16;
end

if nargin<2
    %total order default
total_order=30;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%finding the size of the weight matrix
[n,ncols]=size(d);
%finding number of evaluation points
aiter=length(alphafine);

%%%% block of partial input checking %%%%%%%%%%%%%%%

%making sure alphafine is a column vector

      if aiter==1
       error('alphafine needs to be a column vector') 
    end
    
%transforming a row into a column vector (when necessary)   
      if size(alphafine,2)>1
        alphafine=alphafine';
    end
 
 if length(iter)>1
        error('iter should be a scalar integer')
    end
    
 if length(total_order)>1
       error('total_order should be a scalar integer')
   end
   
   if length(asize)>1
       error('asize should be a scalar')
   end
       
 if n~=ncols
     error('weight matrix needs to be n by n')
   end
   
   if ((randstate<0)|logical( rem(randstate,1)~=0 )|( length(randstate)>1 ) );
       error('randstate must be a non-negative integer scalar');
   end
      
   if max(alphafine)>=1
       error('alphafine should have a maximum element of less than 1 and a minimum element of greater than the inverse of the minimum eigenvalue')
   end
   
   if (asize>=1)|(asize<=0)
       error('asize must lie in (0,1)')
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Measuring traces of 1st two moments. td=[tr(D) tr(D*D)/2]. Note, tr(D) should be equal to 0.
td=full([0;sum(sum(d.*(d')))/2]);

%The scalar total_order specifies total number of moments and is comprised of 2 exact moments, tr(D),tr(D*D)
%and o-2 stochastic moments u'(D^j)u/(u'u).
%This block accumulates stochastic moments from order 1 to total order in columns
%and independent realiztions of the stochastic moments in rows.
mavmomi=zeros(total_order, iter);
for j=1:iter;
u=randn(n,1);
v=u;
utu=u'*u;
for i=1:total_order;
v=d*v;
mavmomi(i,j)=n*((u'*v)/(i*utu));
end;
end;

%this substitutes in the first two exact moments for the first two stochastic moments
mavmomi(1:length(td),:)=td(:,ones(iter,1));

%averages across iterations
avmomi=mean(mavmomi')';

clear u,v; %minor memory savings

%weighting moments by respective powers of a
seq_order=(1:total_order);
oiter=ones(total_order,1)./seq_order';
even_ind=((-1).^seq_order>0)'./seq_order';
srvs=zeros(iter, aiter);
sum_all_weights=zeros(aiter, 1);
sum_even_weights=zeros(aiter,1);
for ii=1:aiter
    wvec=(alphafine(ii).^seq_order);
srvs(:,ii)=-(wvec*mavmomi)'; 
sum_all_weights(ii)=wvec*oiter;
sum_eve_weights(ii)=wvec*even_ind;
end

%Estimated ln|I-aD| using mixture of exact, stochastic moments
%exact from 1 to oexact, stochastic from (oexact+1) to total_order
lndetmat=mean(srvs)'; %average of independent realizations of Krylov space
sderr=(std(srvs)/sqrt(iter))'; %standard error of independent realizations of Krylov space

%lower bound computation
%%The Linear Algebra and its Applications paper uses the following:
%fbound=((n*alpha.^(total_order+1))./((total_order+1)*(1-alpha)))';

%However, for matrices with max(abs(eigenvalue)) equal to or less than 1,
%and an even order_total, the following bound is tighter. 

last_moment=avmomi(end);
posind=(alphafine>=0);

altlowpos=last_moment*(log(1-alphafine)+sum_all_weights);
altlowneg=last_moment*0.5*(log(1-alphafine.*alphafine)+sum_even_weights);
fbound=posind.*altlowpos+(1-posind).*altlowneg;

%confidence limits, with lower limit biased downward (more conservative)
%The paper used a normal approximation or a t for small iter.
%Using Chebyshev's Theorem should yield a very conservative asize% confidence interval. 
cfactor=sqrt(1/asize);

%confidence limits of size asize
low_asize=(lndetmat-cfactor*sderr+fbound);
high_asize=(lndetmat+cfactor*sderr);

%AR parameter, lower confidence limit, estimated log-det, upper confidence limit
confide=[alphafine low_asize lndetmat high_asize];

