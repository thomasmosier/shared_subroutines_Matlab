function [bmax, srds, prhighers, emax, logliks]=fmess_car2(x, y, d)
%
% [bmax, srds, prhighers, emax, logliks]=fmess_car2(x, y, d)
%
%This function computes spatial conditional autoregression when the variance-covariance function is a matrix exponential.
%This specification has many computational gains, and I have computed
%spatial autoregressions using millions of observations with this function.
%
%INPUTS:
%
%x is a n by k matrix of independent variables. 
%
%y is a n by 1 column vector of observations on the dependent variable.
%
%d is a symmetric n by n spatial weight matrix. This function assumes d has a maximum eigenvalue of 1
%
%OUTPUTS:
%
%bmax is the k+1 by 1 vector of parameter estimates with the regression parameter estimates as the first k elements and 
%alphamax, the global spatial dependence parameter (a scalar), as the last element.
%
%srds is the k+1 by 1 vector of signed root deviances. These have a t-ratio like interpretation.
%The first k elements are for the delete-1 hypotheses and the last is for alpha=0. Squaring the srds yields the likelihood ratios.
%
%prhighers is a k+1 by 1 vector of estimates of obtaining higher likelihood
%ratio under repeated sampling.
%
%emax is n by 1 column vector of MESS residuals.
%
%logliks is a iter by (k+1) matrix of profile log-likelihoods.
%
%
%NOTES:
%
%Divide the dependence parameter by 2 to get MESS-SAR in the covariance.
%
%Take the complement (negate) the parameter estimate to get the estimated
%parameter in the inverse variance-covariance matrix.
%
%One can square the signed root deviances to get likelihood ratios, if desired. 
%
%I computed prhigher using likelihood ratios, but one could use the t
%distribution as well.
%
%Theoretically, X'expm(-alpha*D)X is always positive definite. 
%Nevertheless, the Cholesky decomposition can fail with imperfectly conditioned design matrices, especially for large values
%of the dependence parameter alpha. Try rescaling the data (possibly mean-centering -- and avoiding the intercept), 
%and then reducing the limits. Try not to use many similar variables. For example, using square feet of living area,
%total square feet of area, squared-area, number of rooms, number of bedrooms, number of bathrooms, kitchen size, and so forth 
%as independent variables in predicting house prices could lead to an ill-conditioned design matrix.
%You should avoid using data matrices where cond(x'x)>10^10 or so. In addition, this
%is more of a problem when using small n (which usually means k/n is larger). 
%
%If you use the function, please cite:
%
%LeSage, James and R. Kelley Pace, Spatial Dependence in Data Mining,
%Data Mining for Scientific and Engineering Applications, 
%Edited by Robert L. Grossman, Chandrika Kamath, Philip Kegelmeyer, Vipin Kumar, and Raju R. Namburu, 
%Kluwer Academic Publishing, 2001.
%
%Kelley Pace, www.spatial-statistics.com, 12/25/02.

%maximum order of D to use in computing the matrix exponential
%a sufficient condition of weight matrices with maximum eigenvalue of 1 is
%to make (max(abs(alphavec))^(q))/(q)! a small number (eg. 0.005)
%lower this to reduce time, if the max(abs(alphavec)) is sufficiently low.
%I have set this very high to allow the use of larger alpha.
q=24;


%The prob_points gives a grid in terms of the more conventional AR parameter
%You can edit this to improve performance, or to decrease the range which may 
%avoid errors with the Cholesky decomposition.
prob_points=[0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.875 0.9 0.925 0.95 0.96 0.97 0.98 0.99];
%This translates it back to the MESS space
alphaneg=log(1-prob_points);
%This gives both AR and MA parameters -- if you change it make sure alphavec contains a 0 (this yields OLS)
alphavec=[alphaneg 0 -alphaneg];
%This finds the lower and upper bounds of the coarse grid, and given an increment, this determines the fine grid
low_bound=min(alphavec);
up_bound=max(alphavec);
fine_increment=0.001;
alphavec_fine=(low_bound:fine_increment:up_bound);

%number of elements of the coarse grid of spatial dependence parameters
iter=length(alphavec);

%measure the number of observations and the number of variables
[n,k]=size(x);
%augment this by 1 to give the dimension of xymat below
k1=k+1;

%sequence of 1...q
qseq=1:q;
%sequence of inverse factorials (0! ... (q-1)!)
factdivisor=gamma(qseq);

%preallocating storage, and identity matrix
eyer=eye(k);
sseall=zeros(iter,(k+1));
bests=zeros(iter,k);

%place x any y side-by-side
xymat=[x y];   

%This loop finds all moment matrices of [X y]'D^{j}[X y] for j=1...q
dxy=xymat;
for j=1:q
    moms(:,:,j)=dxy'*xymat;
    dxy=d*dxy;
end

%This loop finds the max log-lik for each possible value of alpha in our coarse grid
for i=1:iter;

%we use each possible value of alpha in our grid    
alpha=alphavec(i);

%for a given alpha, finds the vector of powers
alphaparmvec=alpha.^(qseq-1);
%scales these by the sequence of factorials
messw=alphaparmvec./factdivisor;   

%This finds the weighted sum of the individual moment matrices
mommat=zeros(k1,k1);
for j=1:q
    mommat=mommat+moms(:,:,j)*messw(j);
end
 
%We use the Normal Equation approach with Cholesky decompositions.
   xtx=mommat(1:k,1:k);
   cholxtx=chol(xtx);
   xtyi=mommat(1:k,k1);
   eq1=cholxtx'\xtyi;
   b=cholxtx\eq1;
   sse=mommat(k1,k1)-b'*xtx*b;
 
%computes rise in SSE when dropping each variable individually
bz=(cholxtx')\eyer;
bsqterm=((b.^2))./(sum(bz.^2)');

%creates matrix of SSE values - first column unrestricted;next k columns restricted SSE
sseall(i,:)=sse+[0 bsqterm'];

end;
%ends the coarse increment loop with matrix of sse for model plus all delete-1 models over coarse alpha grid
 
%interpolates to a finer grid of values of alpha (alphavec_fine)
sseall_fine=spline(alphavec,sseall',alphavec_fine)';

%OLS sse and loglik
alpha0ind=(alphavec==0);
sseols=sseall(alpha0ind,1);
loglikols=-n*0.5*log(sseols);

%computes MESS-CAR log-likelihood up to a constant
logliks=-(n/2)*log(sseall_fine);

%finds maximum likelihoods and associated indices for each hypothesis
[maxliks,maxindc]=max(logliks);

%finds optimal autoregressive parameters
maxalphas=alphavec_fine(maxindc);

%finds MESS-CAR optimal alpha for the complete model
alphamax=maxalphas(1);

%computes likelihood ratios (twice difference in log-likelihoods) and prob of higher chi2
maxliksall=[maxliks(2:end)';loglikols];%delete-1 models, followed by OLS (alpha=0)
globalmaxlik=max(maxliks);%max lik from overall model
likratios=2*(globalmaxlik-maxliksall);%likelihood ratios
prhighers=1-gammainc(likratios/2,1/2);

%This takes the MESS-SAR estimated alpha (which equals 0.5*MESS-CAR estimated alpha)
%and finds the appropriate powers of alpha
maxalphaparmvec=(0.5*alphamax).^(qseq-1);

%This reweights the powers with the inverse of the factorial
maxmessw=maxalphaparmvec./factdivisor;   

%This loop computes expm(-alpha*0.5*D)[X y]
dxy=xymat;
mxymat=zeros(n,k1);
for j=2:(q+1)  
    mxymat=maxmessw(j-1)*dxy+mxymat;
    dxy=d*dxy;
end

%giving expm((-alpha*0.5)*D[X y] identifying names for transformed X,y
xmax=mxymat(:,1:k);
ymax=mxymat(:,end);

%finding MESS-CAR==MESS-SAR parameter estimates
bmax=xmax\ymax;

%This finds the residuals using the ML estimate in the transformed space
%This should be equivalent to MESS-CAR residuals as well.
%One can test emax, since emax'*emax should approximately equal exp(-(2/n)*max(max(loglik)))
emax=(ymax-xmax*bmax);

%Stacking estimates. Also, this returns the MESS parameter in the variance-covariance matrix, as
%opposed to the inverse variance-covariance matrix (by taking complement). 
bmax=[bmax;-alphamax];

%For interpretability, and scaling we convert the likelihood ratio into a signed-root-deviance.
%This has a t-ratio like interpretation
srds=sqrt(likratios).*sign(bmax);
