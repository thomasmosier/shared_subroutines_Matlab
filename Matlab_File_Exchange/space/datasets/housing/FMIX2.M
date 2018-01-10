function [bmax, srds, prhighers, emax, logliks]=fmix2(x, y, d, detvalz)
%
%  [bmax, srds, prhighers, emax, logliks]=fmix2(x, y, d, detvalz)
%
%This function computes either autoregressive or mixed regressive spatially autoregressive estimates. 
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
%
%OUTPUTS:
%
%bmax is the k+1 by 1 vector of parameter estimates with the regression parameter estimates as the first k elements and 
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
%Pace, R. Kelley, and Ronald Barry, “Quick Computation of Regressions with a Spatially Autoregressive Dependent Variable,” 
%Geographical Analysis, Volume 29, Number 3, July 1997, p. 232-247.
%
%Kelley Pace, www.spatial-statistics.com, written on 6/23/97, revised on 12/25/02.


%supplies Chebyshev approximation when user does not supply log-determinants
if nargin<4
    [detvalz]=fsubdet_chebyshev2(d);
end

[n,k]=size(x);
[aiter,dcols]=size(detvalz);

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
    error('This appears to be neither a log-determinant matrix (aiter by 2) or bounded log-determinant matrix (aiter by 4).')
end

%defining some useful variables for later uses
o=ones(aiter,1);
eyer=eye(k);
alphafinemat=alphafine(:,ones(1,k));

%the following compute components of the moment matrices
dy=d*y;
xtx=x'*x;
xty=x'*y;
xtdy=x'*dy;
yty=y'*y;
ytdy=y'*dy;
ytd2y=dy'*dy;

%computing component parameter estimates
 cholxtx=chol(xtx);
 bz=(cholxtx')\eyer;
 invbz2=1./(sum(bz.^2));
 
 part1=cholxtx'\xty;
 b=cholxtx\part1;
  
 part2=cholxtx'\xtdy;
 bd=cholxtx\part2;
 
%forming components of global sse
 sseyy=yty-b'*xtx*b; 
 ssedydy=ytd2y-bd'*xtx*bd;
 sseydy=ytdy-b'*xtx*bd;
 
 sseglobal=sseyy(o,:)-2*alphafine.*sseydy(o,:)+(alphafine.*alphafine.*ssedydy(o,:));
 sseglobalexpand=sseglobal(:,ones(1,k));
  
%mixed regressive spatially autoregressive estimate 
 bdalpha=(b(:,o')'-alphafine*bd');
 bdalpha2=bdalpha.^2;
 
 %tricks using restricted least-squares to compute the rise in SSE from delete-1 restrictions
 invbz2expand=invbz2(o,:);
 delsse=invbz2expand.*bdalpha2;

%unresticted and restricted SSEs
ssefines=[sseglobal (sseglobalexpand+delsse)];

%finding optima
[amaxs, maxliks, logliks, maxindc]=famax(logdets, ssefines, alphafine, n);
alphamax=amaxs(1);

%using optimum global estimates to compute residuals, global parameters
ymax=(y-alphamax*dy);
part1max=cholxtx'\(xty-alphamax*xtdy);
bmax=cholxtx\part1max;
emax=ymax-x*bmax;

%here is OLS logliklihood for comparison
loglikols=(-n*0.5)*log(yty-b'*xtx*b);

%inference exact or likelihood dominance
if exact_inference==1;
%computes likelihood ratios for beta and alpha
likratios=(2*abs(max(maxliks)-maxliks(2:end)))';
likratioalpha=2*(maxliks(1)-loglikols);
else
%likelihood dominance route
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
%prob of higher test statistic
prhighers=1-gammainc(totlikratios/2,1/2);
%signed root deviance
srds=sqrt(totlikratios).*sign(bmax);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [amaxs, maxliks, logliks, maxindc]=famax(logdetfine, ssefines, alphafine, n)
[iter,k1]=size(ssefines);
ok1=ones(1,k1);
logliks=logdetfine(:,ok1)-0.5*n*log(ssefines);%SAR likelihood function
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

