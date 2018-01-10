function [bmax, srds, prhigher, emax, maxlik]=fmess_ar2(x, y, d)
%
% [bmax, srds, prhigher, emax, maxlik]=fmess_ar2(x, y, d)
%
%This function estimates an autoregressive matrix exponential spatial
%specification (MESS).
%
%
%INPUT:
%
%x contains n observations on k independent variables (n by k)
%
%y contains n observations on a single dependent variable (n by 1)
%
%d is a n by n spatial weight matrix (n by n)
%
%
%OUTPUT:
%
%bmax is the k+1 by 1 vector of parameter estimates with the SAR regression parameter estimates as the first k elements and 
%alphamax, the global spatial dependence parameter (a scalar), as the last element.
%
%srds is the k+1 by 1 vector of signed root deviances associated with bmax. These have a t-ratio like interpretation. 
%The first k elements are for the delete-1 hypotheses and the last is for alpha=0. Squaring the srds yields the likelihood ratios.
%
%prhighers is a k+1 by 1 vector of estimates of obtaining a higher likelihood ratio under repeated sampling.
%
%emax is n by 1 column vector of MESS residuals.
%
%maxlik is a scalar giving the maximum of the profile log-likelihood
%
%
%NOTES:
%
%If you have extremely large spatial dependence parameters such as alphamax
%of over 6, you may need to increase nq.
%
%If you use this function, please cite: 
%
%You can see the following for more information:
%
%LeSage, James and R. Kelley Pace, “Spatial Dependence in Data Mining,” 
%Data Mining for Scientific and Engineering Applications, 
%Edited by Robert L. Grossman, Chandrika Kamath, Philip Kegelmeyer, Vipin Kumar, and Raju R. Namburu, 
%Kluwer Academic Publishing, 2001.
%
%Written by Kelley Pace, www.spatial-statistics.com, 12/24/02.

%getting the number of observations, variables
[n,k]=size(x);
%forming an identity matrix for later use
eyer=eye(k);
%initializing restricted SSEs for later use
qrest=zeros(k,1);
%highest term (less 1) in matrix exponential Taylor series
nq=12;

%this block creates the Krylov subspace associated with the weight matrix d
wiy=y;
ymat=zeros(n,nq);
for i=1:nq;
ymat(:,i)=wiy;
wiy=d*wiy;
end;

%form nq by nq matrix of sum of squares
ssqmat=ymat'*ymat;
%form k by k moment matrix
xtx=x'*x;
%form k by nq moments
xty=x'*ymat;

%Cholesky decomposition of moment matrix
cholxtx=chol(xtx);
%The next two steps solve for the least squares estimates
eq1=cholxtx'\xty;
bmat=cholxtx\eq1;

%This forms the overall nq by nq matrix containing the sum of squared errors
ssemat=ssqmat-bmat'*xtx*bmat;

%This calls the closed-form MESS estimator which returns a vector containg
%weights given to the Krylov space and the optimal sum-of-squared errors.
[wpstarvec, qstar]=fmessroot(ssemat);

%This function naturally estimates exp(-0.5*D). So we multiply by -2.0 to
%obtain the estimated covariance function dependence parameter
 alphamax=-2.0*wpstarvec(2);   

%bmax is the MESS maximum likelihood estimate of beta augmented with alphamax
 bmax=bmat*wpstarvec;
 
 %These are the MESS residuals
 emax=(ymat*wpstarvec)-x*bmax;

 %These terms are used in computing the restricted sum-of-squared errors
bz=(cholxtx')\eyer;
restinv=1./sum(bz.^2);

%Each loop computes the SSE matrix for each delete-1 submodel, and optimal
%MESS restricted SSE
for jhyp=1:k
b_jhyp=bmat(jhyp,:);
sse_rest=ssemat+(b_jhyp'*b_jhyp)*restinv(jhyp);
[prestvec, qrest(jhyp)]=fmessroot(sse_rest);
end

%here is OLS logliklihood for comparison
loglikols=(-n*0.5)*log(ssemat(1,1));

%This compute the overall MESS maximum likelihood statistic
maxlik=-0.5*n*log(qstar);

%This computes the delete-1 restricted MESS maximum likelihood statistics,
%and alphamax=0.
logliks=[-0.5*n*log(qrest);loglikols];

%This forms the deviances
likratios=2*(maxlik-logliks);

%This gives the probability of observing a higher likelihood ratio under
%repeated sampling
prhigher=1-gammainc(likratios/2,1/2);

%The signed root deviance is the square root of the deviance given the sign
%of the parameter
bmax=[bmax;alphamax];
srds=sqrt(likratios).*sign(bmax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wpstarvec, qstar]=fmessroot(ssemat)
%
%This computes the MESS closed-form solution
%
%INPUT:
%
%SSEMAT is a nq by nq matrix containing the SSEs 
%
%OUTPUT:
%
%wpstarvec is the optimal weights to give to the columns of ymat, bmat.
%
%qstar is a scalar optimum overall SSE for a model
%
%Kelley Pace, 12/24/02

%Getting dimension information and forming matrix exponential Taylor series
%weights
nq=size(ssemat,1);
nq1=nq-1;
rrr=(1./gamma(1:nq))';
drrr=diag(rrr);

%Applying the deterministic weights to the SSE matrix
zzz=drrr*ssemat*drrr;

%We manipulate the matrix to isolate all polynomials of a given degree
flipzzz=flipud(zzz);

%We sum across all same degree terms
   p=fliplr(sum(spdiags(flipzzz,-nq1:nq1)));
   
   %We find the roots
   posstar=roots(polyder(p));
   
   %We isolate the real root, which should be optimal
   realind=imag(posstar)==0;
   pstar=posstar(realind);

   %We reform the MESS weights given the optimal parameter
   pstarvec=(pstar.^(0:nq1))';
   
   %This allows computation of the optimal overall SSE
   qstar=pstarvec'*zzz*pstarvec;
   
   %This gives the optimal weighting of the columns of bmat, ymat
   wpstarvec=drrr*pstarvec;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 