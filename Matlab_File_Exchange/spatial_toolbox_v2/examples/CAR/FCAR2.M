function [bmax, srds, prhighers, emax, logliks]=fcar2(x, y, d, detvalz, itersub)
%
%  [bmax, srds, prhighers, emax, logliks]=fcar2c(x, y, d, detvalz, itersub)
%
%This function computes conditional spatial autoregression (CAR or CSG). 
%
%
%INPUTS:
%
%x is a n by k matrix of independent variables. 
%
%y is a n by 1 column vector of observations on the dependent variable.
%
%d is a symmetric n by n spatial weight matrix. This function assumes d has a maximum eigenvalue of 1.
%
%detvalz is either a iter by 2 matrix containing a vector of autoregressive parameters in the first column and 
%the associated log-determinant in the second column or a iter by 4 matrix containing a vector of autoregressive
%parameters in the first column, the associated lower bound on the log-determinant in the second column, 
%the associated log-determinant in the third column, and the associated lower bound on the log-determinant in the fourth column.
%
%itersub is the number of iterations used in interpolating the sse function. As this becomes larger the results approach the exact
%ones.
%
%
%OUTPUTS:
%
%bmax is the k+1 by 1 vector of parameter estimates with the CAR regression parameter estimates as the first k elements and 
%alphamax, the global spatial dependence parameter (a scalar), as the last element.
%
%srds is the k+1 by 1 vector of signed root deviances. These have a t-ratio like interpretation. If the exact log-determinant 
%(or a very close approximation) is used, this has a standard interpretation. For less precise approximations, the srds is 
%based on likelihood dominance inference and should provide lower bounds to the exact srds.
%
%prhighers is a k+1 by 1 vector of estimates of obtaining a higher likelihood ratio under repeated sampling.
%
%emax is n by 1 column vector of residuals in the spatially transformed space.
%
%loglik is a iter by (k+2) matrix of profile log-likelihoods (first column full model, remaining columns delete-1 submodels).
%
%
%NOTES:
%
%One can square the signed root deviances to get likelihood ratios, if desired. 
%
%This function uses some of the techniques discussed in:
%
%Pace, R. Kelley, and James LeSage, “Likelihood Dominance Spatial Inference,” forthcoming in Geographical Analysis.
%
%Pace, R. Kelley, and Ronald Barry, “Fast CARs,” 
%Journal of Statistical Computation and Simulation, Volume 59, Number 2, 1997, p. 123-147.
%
%Kelley Pace, www.spatial-statistics.com, written on 6/23/97, revised on 12/25/02.



if nargin<5
%itersub gives how many interpolation points are used
itersub=20;
end

%supplies Chebyshev approximation when user does not supply log-determinants
if nargin<4
    [detvalz]=fsubdet_chebyshev2(d);
end

[n,k]=size(x);
[aiter,dcols]=size(detvalz);
%NOTE: actually aiter is iter with respect to the description of the variables in the header. 

%this gives names to the autoregressive parameter vector and the log-determinant vector
alphafine=detvalz(:,1);
if dcols==2
logdets=detvalz(:,2);
exact_inference=1;
elseif dcols==4
    logdets=detvalz(:,3);
    lowerbounds=detvalz(:,2);
    upperbounds=detvalz(:,4);
    exact_inference=0;
else
    error('This appears to be neither a log-determinant matrix (iter by 2) or bounded log-determinant matrix (iter by 4).')
end

%the following compute components of the moment matrices
xtx=x'*x;
xty=x'*y;

dy=d*y;
dx=d*x;

xtdx=x'*dx;
xtdy=x'*dy;

yty=y'*y;
ytdy=y'*dy;

%finding limits of input dependence parameter
low_bound=min(alphafine);
up_bound=max(alphafine);

if low_bound>=0
%The log-determinant is nonlinear, and so this spaces evaluation points equally for the function (appproximate).
alphacoarse=1-exp(linspace(0,log(0.01),(itersub+1)));
else
alphacoarsepos=1-exp(linspace(0,log(0.01),(itersub+1)));
alphacoarseneg=alphacoarsepos*low_bound;
alphacoarse=[alphacoarseneg alphacoarsepos(2:end)];
end

coarseiter=length(alphacoarse);
eyer=eye(k);
sseall=zeros(coarseiter,(k+1));

for i=1:coarseiter;
   
   alpha=alphacoarse(i);
   
   %computes moment matrices from their components
   xtxi=(xtx-alpha*xtdx);
   xtyi=(xty-alpha*xtdy);
   ytyi=(yty-alpha*ytdy);
   
   %takes cholesky decomposition
   cholxtx=chol(xtxi);
   eq1=cholxtx'\xtyi;
   b=cholxtx\eq1;
   sse=ytyi-b'*xtxi*b; 
   
   %computes rise in SSE when dropping each variable individually
bz=(cholxtx')\eyer;
bsqterm=((b.^2))./(sum(bz.^2)');

%creates matrix of SSE values - first column unrestricted;next k columns restricted SSE
sseall(i,:)=sse+[0 bsqterm'];

end;
%ends the coarse increment loop with matrix of sse for model plus all delete-1 models over coarse alpha grid

%interpolates to a finer grid of values of alpha (alphafine)
ssefines=spline(alphacoarse,sseall',alphafine')';

%obtaining estimates
[amaxs, maxliks, logliks, maxindc]=famax(logdets, ssefines, alphafine, n);
alphamax=amaxs(1);
ymax=(y-alphamax*dy);
xmax=(x-alphamax*dx);
bmax=(xtx-alphamax*xtdx)\(xty-alphamax*xtdy);
emax=ymax-xmax*bmax;

%here is OLS logliklihood for comparison
bols=xtx\xty;
loglikols=(-n*0.5)*log(yty-bols'*xtx*bols);

%inference exact or likelihood dominance
if exact_inference==1;
%computes likelihood ratios for beta and alpha
likratios=(2*abs(max(maxliks)-maxliks(2:end)))';
likratioalpha=2*(maxliks(1)-loglikols);
else
%likelihood dominance inference
ssefull=ssefines(:,1);
%differences in sse from global to each delete-1 restriction
restdiffs=-n*(log(ssefull(:,ones(k,1)))-log(ssefines(:,2:end)));

%obtaining bounds on alpha for each delete-1 restriction
[amaxslow, maxlikslow, loglikslow, maxindclow]=famax(lowerbounds, ssefines, alphafine, n);
[amaxsup, maxliksup, logliksup, maxindcup]=famax(upperbounds, ssefines, alphafine, n);

%loop to compute minimum difference in SSE over restricted range for each delete-1 restriction
lbind=zeros(k,1);
ubind=zeros(k,1);
likratios=zeros(k,1);
for iii=1:k
lbind(iii)=maxindclow(iii+1);
ubind(iii)=maxindcup(iii+1);
restdiffssub=restdiffs(lbind(iii):ubind(iii),iii);
[likratios(iii), rinds]=min(restdiffssub);
end

%finds minimum likelihood ratio for alpha based on restrictions
likratioalpha=2*(min(maxlikslow(1),maxliksup(1))-loglikols);
end

%putting estimates together
bmax=[bmax;alphamax];
totlikratios=[likratios;likratioalpha];
prhighers=1-gammainc(totlikratios/2,1/2);
srds=sqrt(totlikratios).*sign(bmax);%SRD is signed root deviance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [amaxs, maxliks, logliks, maxindc]=famax(logdetfine, ssefines, alphafine, n)
[iter,k1]=size(ssefines);
ok1=ones(1,k1);
logliks=0.5*logdetfine(:,ok1)-0.5*n*log(ssefines);%CAR likelihood function
[maxliks, maxindc]=max(logliks);
amaxs=alphafine(maxindc);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [detbounds]=fsubdet_chebyshev2(d, alphafine)
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

