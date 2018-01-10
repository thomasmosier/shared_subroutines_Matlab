function [eouts, erecs, alphamaxes, bmaxes]=fwholesale2(xall, yall, d,  start, nsub, obsindices, xcoord, ycoord)
%[eouts,erecs,alphamaxes,bmaxes]=fwholesale2(xall, yall, d,  start, nsub, obsindices, xcoord, ycoord)
%
%This function computes a set of spatial autoregressive local
%estimates (SALE) where the ith SALE are centered at xcoord(obsindices(i)),
%ycoord(obsindices(i)). Each SALE uses at minimum start observations and at
%maximum nsub observations. 
%
%INPUT:
%
%xall is a n by k matrix containing the n observations on k independent variables for the
%global sample
%
%yall contains y and its spatial lag dy where dy=d*y . This is a n by 2
%matrix.
%
%d is an n by n symmetric spatial weight matrix
%
%start is the minimum number of observations used for any SALE (start>k)
%
%nsub is the maximum number of observations used for any SALE
%
%obsindices are the observations numbers of the desired centers of the SALE
%samples. This is a jiter by 1 vector. 
%
%xcoord, ycoord are both n by 1 vectors of locational coordinates
%
%
%OUTPUT:
%eouts are nsub residuals of the holdout observation at
%xcoord(obsindices(i)), ycoord(obsindices(i)) for i=1....jiter. It is a
%nsub by jiter matrix.
%
%erces are nsub residuals from predicting the i+1 observations based upon observations 2 to i (observation 1 is held out as well).
%These are recursive residuals
%
%alphamaxes are the matrix of nsub by jiter spatial dependence parameters.
%
%bmaxes are a sequence of beta estimates for the sequence of samples
%start:i for i=start:nsub.
%
%NOTES:
%
%If you wish to run SALE for a large number of indices, you may need to
%edit this function to not store every statistic for each subsample
%regression.
%
%If you use these functions, please cite:
%Pace, R. Kelley, and James LeSage, “Spatial Autoregressive Local Estimation,” 
%Spatial Statistics and Spatial Econometrics, Edited by Art Getis, Palgrave, 2003.
%
%Written by Kelley Pace, www.spatial-statistics.com, 12/21/02

%We begin by obtaining the Chebyshev coefficients and the spatial
%dependence grid.
[comboterm, alphavec]=fchebyshev_initialize2(d);

%implied number of iterations
jiter=length(obsindices);

%size of the design matrix
[n,k]=size(xall);

%independent and dependent variables
alldat=[xall yall];

%preallocating
erecs=zeros(nsub,jiter);
eout=zeros(nsub,jiter);
maxliks=zeros(jiter,1);
alphamaxes=zeros(nsub,jiter);
bmaxes=zeros(nsub,k,jiter);

for i=1:jiter
    
%finding observation associated with index i   
j=obsindices(i);

%computes squared Eucliean distance away from initial point
d2=sqrt((xcoord-xcoord(j)).^2+(ycoord-ycoord(j)).^2);
%sorts, and finds associated indices -- the squaring of the distance
%doesn't affect the ordering
[d2s,dind]=sort(d2);

%holdout at initial observation 1
dindout=dind(1);
%subsample of initial observation
subdatout=alldat(dindout,:);
%taking subsample and partitioning it into the independent variables, dependent variable, and lagged dependent variable
xout=subdatout(1,1:k);
yout1=subdatout(1,(k+1));
dyout1=subdatout(1,(k+2));

%takes nsub observations 2 to (nsub+1) ordered by Euclidean distance away from observation 1
dindsub=dind(2:(nsub+1));
%subsample extending from initial observation
subdat=alldat(dindsub,:);
%taking subsample and partitioning it into the independent and dependent variable and lagged dependent variable
xt=subdat(:,1:k);
yt=subdat(:,(k+1):(k+2));

%local weight matrix -- global with rows, columns pertaining to subsample separated out
%This weight matrix has neighbors outside the subsample. Thus, the spatial parts of yt and xt 
%may extend outside the subsample. However, only simultaneity in the sample matters.
wswj=d(dindsub,dindsub);

%sequence of log-determinants using local weight matrix wswj
%[detvalzt,alphavec]=fdet_updates2(wswj);
%or use the much faster, approximate solution
%[lowerbounds, detvalzt, upperbounds, alphavec]=fchebyshev_updates2c(wswj);
[detvalzt]=fchebyshev_updates2(wswj, comboterm);

%gets individual SALE parameter estimates and recursive residuals
[alphamax,bmax,erecurspace]=fmixsale2(xt,yt,detvalzt,alphavec,start);

%saving sequence of optimal alpha, optimal beta
alphamaxes(:,i)=alphamax;
bmaxes(:,:,i)=bmax;

% saving sequence of recursive residuals
erecs(:,i)=erecurspace;

%calculates sequence of residuals when predicting observation 1 (for each local estimation attempt)
eouts(:,i)=yout1-alphamax*dyout1-bmax*xout';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [comboterm, rhovec]=fchebyshev_initialize2(d)
% [comboterm, rhovec]=fchebyshev_initialize2(d)
%
%This function provides the Chebyshev coefficients for fchebyshev_updates2 which computes the 
%log-determinant of (I-aD(1:i,1:i)) for i=1...n, where a can go from -1 to 1 (for this routine). The routine
%assumes the input matrix d is symmetric with a trace of 0.
%
%INPUT:
%d, a n by n symmetric weight matrix with zeros on the diagonal
%
%OUTPUT:
%comboterm is an iter by 3 matrix comprised of the Chebyshev constant,
%linear, and quadratic coefficients.
%
%rhovec, an iter by 1 vector comprised of the values of the spatial
%dependence parameter.
%
%
%If you use this, please cite:
%
%Pace, R. Kelley, and James LeSage, 
%"Chebyshev Approximation of Log-determinants of Spatial Weight Matrices,"
%Computational Statistics and Data Analysis, forthcoming.
%
%
%Kelley Pace, December 7, 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%we need the size of the matrix for the approximation
[nn,nn]=size(d);

%we need a row unit vector
o=ones(1,nn);

%we need a row sequence vector
nnseq=(1:nn);

%td1=trace(d); If one needs to compute this. Usually is zero by
%construction. The Taylor bounds assume this.
td1=0;

%trace(D(1:i)*D(1:i)) for sequence i=1 ... n
d_upper_SSQ=2.0*sum(triu(d).^2);
td2_sequence=cumsum(d_upper_SSQ);

%spatial dependence parameters -- edit these to taste
%rhovec=((-99:99)/100)';%For both positive and negative dependence
rhovec=((0:99)/100)'; %For positive dependence.

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

for i=1:length(rhovec); 
rho=rhovec(i);   
for j=1:nposs;
cposs(i,j)=(2/nposs)*sum(log(1-rho*xk).*cos(pi*(j-1).*(seq1nposs-0.5)/nposs));
end;
end;

%This multiplies the Chebyshev coefficients times the polynomials
comboterm=cposs*rezzy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chebyshevest]=fchebyshev_updates2(d, comboterm)
% [chebyshevest]=fchebyshev_updates2(d, comboterm)
%
%This function returns a sequence of quadratic lower Taylor series bounds, quadratic Chebyshev approximations, and upper
%Taylor series bounds for the log-determinant of (I-aD(1:i,1:i)) for i=1...n, where a can go from -1 to 1 (for this routine). The routine
%assumes the input matrix d is symmetric with a trace of 0.
%
%INPUT:
%d, a n by n symmetric weight matrix with zeros on the diagonal
%comboterm, an iter by 3 matrix of constant, linear, and quadratic
%coefficients needed for the Chebyshev approximation.
%
%OUTPUT:
%chebyshevest is a iter by 1 matrix comprised of a sequence of Chebyshev
%log-determinant approximation (ln|I-aD(1:i,1:i)| for i=1...n.
%
%
%If you use this, please cite:
%
%Pace, R. Kelley, and James LeSage, 
%"Chebyshev Approximation of Log-determinants of Spatial Weight Matrices,"
%Computational Statistics and Data Analysis, forthcoming.
%
%Kelley Pace, December 7, 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%we need the size of the matrix for the approximation
[nn,nn]=size(d);

%we need a row unit vector
o=ones(1,nn);

%we need a row sequence vector
nnseq=(1:nn);

%td1=trace(d); If one needs to compute this. Usually is zero by
%construction. The Taylor bounds assume this.
td1=0;

%trace(D(1:i)*D(1:i)) for sequence i=1 ... n
d_upper_SSQ=2.0*sum(triu(d).^2);
td2_sequence=cumsum(d_upper_SSQ);

%This is a sequence of size, trace(D), trace(D*D)-n/2] 
%The last vector is also adjusted by -0.5*nnseq which appears in the approximation formula.
lastvec=td2_sequence-0.5*nnseq;
tdvec=[nnseq;td1*o;lastvec];

%This computes the actual Chebyshev quadratic approximation of the
%log-determinant for each value of rho in the grid and for the sequence
%1...nn. lnest should have size iter by nn
chebyshevest=comboterm*tdvec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alphamax,bmax,erecurspace]=fmixsale2(xt,yt,detvalzt,alphavec,start);
%
%[alphamax,bmax,erecurspace]=fmixsale2(xt,yt,detvalzt,alphavec,start)
%
%PURPOSE:
%
%This function calculates SALE (Spatial Autoregressive Local Estimates) given a set of independent and dependent
%variables, a sequence of log-determinants, the grid of spatial dependence
%parameters associated with the log-determinant, and a minimum number of
%observations to be used.
%
%INPUT:
%
%xt represents nsub observations on k independent variables
%
%yt represents nsub observations on the dependent variable
%
%detvalzt represents nsub (or iter) log-determinants for ndets values of
%the spatial dependence parameter found in alphavec
%
%alphavec is a ndet by 1 vector containing the values of the spatial
%dependence parameter alpha
%
%start is the minimum number of observations to use for SALE and this must
%exceed k.
%
%
%If you use these functions, please cite:
%Pace, R. Kelley, and James LeSage, “Spatial Autoregressive Local Estimation,” 
%Spatial Statistics and Spatial Econometrics, Edited by Art Getis, Palgrave, 2003.
%
%Written by Kelley Pace, 12/21/02


%calls to get the recursive errors and regression parameter estimates
[erecur,sset,brecur1,brecur2]=fsseseq3(xt,yt,start);

%takes square of values in the grid of alphas
alphavec2=alphavec.^2;
%measures number of log-determinants
ndets=length(alphavec);
%a constant vector ndets by 1 
odets=ones(ndets,1);
%also ndets by iter after transpose
ssemat=sset';

%finding iter and k
[iter,k]=size(brecur1);

%preallocating
alphamax=zeros(iter,1);
bmax=zeros(iter,k);
loglikmax=zeros(iter,1);
ssemax=zeros(iter,1);
erecurspace=zeros(iter,1);

%finds moments for each a
for i=start:iter
    e2o=ssemat(1,i);
    edo=ssemat(2,i);
    e2d=ssemat(3,i);
    lndet=detvalzt(:,i);
    
%finds log of sequence of variances
logsse=log((odets*e2o-2*alphavec*edo+alphavec2*e2d)/i);

%finds sequence of logliks
loglik=(lndet-(i/2)*logsse);

%finds max-lik
[mxloglik,logind]=max(loglik);

%maximum likelihood estimate of alpha
alphamax(i)=alphavec(logind);

%recursive maximum likelihood estimate of regression parameters
bmax(i,:)=brecur1(i,:)-alphamax(i)*brecur2(i,:);
%recursive variance at maximum
ssemax(i)=(e2o-2*alphamax(i)*edo+alphamax(i)*alphamax(i)*e2d)/i;
%recursive AR residuals
erecurspace(i)=erecur(i,1)-alphamax(i)*erecur(i,2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [erecur,sset,brecur1,brecur2]=fsseseq3(xt,yt,start);
%
%[erecur,sset,brecur1,brecur2]=fsseseq3(xt,yt,start)
%
%PURPOSE:
%
%This function estimates recursive residuals, sum-of-squared errors,
%regression parameters.
%
%INPUT:
%
%xt represents nsub observations on k independent variables
%
%yt represents nsub observations on the dependent variable and its spatial
%lag. Hence, yt is a nsub by 2 matrix.
%
%start gives the minimum number of observations to use before beginning
%recursion
%
%OUTPUT:
%
%erecur represents nsub recursive residuals from the prediction of yt and
%its spatial lag. Hence, erecur is a nsub by 2 matrix.
%
%sset is a sequence of sum-of-squared residuals from y, dy, and the
%cross-products of the residuals from predicting y, dy.
%
%brecur1 is the nsub by k recursive regression parameter estimates from
%estimating using y.
%
%brecur1 is the nsub by k recursive regression parameter estimates from
%estimating using dy (spatially lagged y).
%
%
%If you use these functions, please cite:
%Pace, R. Kelley, and James LeSage, “Spatial Autoregressive Local Estimation,” 
%Spatial Statistics and Spatial Econometrics, Edited by Art Getis, Palgrave, 2003.
%
%Written by Kelley Pace, 12/21/02

[n,k]=size(xt);
xtsub=xt(1:(start-1),:);
ytsub=yt(1:(start-1),:);
momtsub=xtsub'*xtsub;
invmomtsub=inv(momtsub);
btsub=invmomtsub*(xtsub'*ytsub);
etsub=ytsub-xtsub*btsub;

bj=btsub;
sj=invmomtsub;
momt=momtsub;
y2t=ytsub'*ytsub;

brecur1=zeros(n,k);
brecur2=zeros(n,k);
dj=zeros(n,1);
erecur=zeros(n,2);
sset=zeros(n,3);

for i=start:n;

bj1=bj;
sj1=sj;

xj=xt(i,:)';
yj=yt(i,:)';
prej=xj'*bj1;
ej=yj'-prej;
erecur(i,:)=ej;

momt=momt+xj*xj';
y2t=y2t+yj*yj';

part1=sj1*xj;
part2=(1+part1'*xj);

dj(i)=sqrt(part2);
bj=bj1+(part1*ej)/part2;
sj=sj1-(part1*part1')/part2;
brecur1(i,:)=bj(:,1)';
brecur2(i,:)=bj(:,2)';
ssetmat=y2t-bj'*momt*bj;
sset(i,:)=[ssetmat(1,1) ssetmat(1,2) ssetmat(2,2)];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [detvalz, alphavec]=fdet_updates2(d);
%
%This function computes a matrix of log-determinants associated with a vector of AR parameters (alpha)
%
%INPUT:
%
%d is a n by n symmetric spatial weight matrix
%
%OUTPUT:
%
%detvalz is an iter by n matrix containing a sequence of log-determinants 
%
%alphavec is a iter by 1 vector of associated spatial dependence parameters
%
%If you use these functions, please cite:
%Pace, R. Kelley, and James LeSage, “Spatial Autoregressive Local Estimation,” 
%Spatial Statistics and Spatial Econometrics, Edited by Art Getis, Palgrave, 2003.
%
%Written by Kelley Pace, 12/21/02

%these settings help optimize the sparse matrix routines
spparms('tight');
spparms('autommd',0);

[n,n]=size(d);
s1=speye(n);

%selecting points to use for the interpolation
iter=100;
alphavec=(((1:iter)-1)/iter)';

detvalz=zeros(iter,n);
for i=1:iter;
alpha=alphavec(i);
u=chol(s1-alpha*d);
%cholesky decomposition
detvalz(i,:)=2*cumsum(log(diag(u)))';
%the cholesky takes the pivots to the 1/2 power. Hence we multiply by 
%2 after taking the sum of the log of the pivots 

end;




                                                                                    