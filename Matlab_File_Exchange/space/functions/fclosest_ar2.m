function [bmax, srds, prhigher,emax, maxlik]=fclosest_ar2(x, y, clo)
%
%[alphamax, loglik, emax, bmax, srds, prhigher]=fclosest_ar2(x, y, closest)
%
%INPUTS:
%
%The matrix x represents n observations on k independent variables.
%
%The vector y represents n observations on the dependent variable.
%
%closest represents either a list of closest neighbor indices such as from fclosest_nn2 
%or a sparse adjacency matrix (s1 when loading smats.mat from fneighbors2).
%
%OUTPUT:
%
%alphamax is a scalar measuring spatial dependence 
%
%maxlik is a scalar giving the maximum of the profile log-likelihood
%
%emax is a n by 1 vector of closest neighbor residuals
%
%bmax is a k by 1 vector of regression parameter estimates
%
%srds is a k by 1 vector of the signed root deviances associated with bmax
%
%prhigher is the probability in repeated sampling of observing a higher srd
%
%NOTES:
%
%Closest neighbor dependence was studied in:
%
%Pace, R. Kelley, and Dongya Zou, 
%“Closed-Form Maximum Likelihood Estimates of Nearest Neighbor Spatial Dependence,” 
%Geographical Analysis, Volume 32, Number 2, April 2000, p. 154-172.
%
%Written by Kelley Pace, www.spatial-statistics.com, 1999, revised
%12/26/02.

%number of observations of dependent variable y
[n, k]=size(x);

[n1, n2]=size(clo);

if (n1==n)&(n2==n)
    %This means we have an n by n weight matrix
    s1=clo;
    %finds number of symmetric pairs used in computing determinant
    %(tr(s1*s1)/2
    c=sum(sum(s1'.*s1))*0.5;
    %computes spatially lagged y
    dy=s1*y;
elseif (n1*n2==n)
    %This means we have a n element vector of indices
nnlist=clo;
%finds number of symmetric pairs used in computing determinant
%if the ith observation has j as neighbor and the jth observation has i
%as a neighbor, we will get i back. Hence, this is a symmetric observation
%by testing which observations have i in the ith position, we can determine
%the number of symmetric elements which when divided by two gives the number
%of symmetric pairs
c=sum(nnlist(nnlist)==(1:n)')/2;

%computes spatially lagged y
dy=y(nnlist);
else
    disp('Error -- neither nnlist or s1 supplied');
end

%puts y and spatially lagged y in a matrix
yall=[y dy];

% moment matrix
xtx=x'*x;

%cholesky decomposition of moment matrix
cholxtx=chol(xtx);
   
%x-y moments
xty=x'*yall;

%computes regression with y and spatially lagged y simultaneously
ball=cholxtx\(cholxtx'\xty);

%finds n by 2 matrix of errors 
eall=yall-x*ball;

%computes 4 by 4 matrix of sum-of-squared errors
ssemat=eall'*eall;

[maxlik, alphamax]=closeopt(ssemat, c, n);

%computing ML  estimate
bmax=ball(:,1)-alphamax*ball(:,2);

%computing associated error
emax=eall(:,1)-alphamax*eall(:,2);

%These terms are used in computing the restricted sum-of-squared errors
eyer=eye(k);
bz=(cholxtx')\eyer;
restinv=1./sum(bz.^2);

%Each loop computes the SSE matrix for each delete-1 submodel, and optimal
%MESS restricted SSE
closeloglik=zeros(k,1);
closealpha=zeros(k,1);
for jhyp=1:k
b_jhyp=ball(jhyp,:);
sse_rest=ssemat+(b_jhyp'*b_jhyp)*restinv(jhyp);
[closeloglik(jhyp), closealpha(jhyp)]=closeopt(sse_rest, c, n);
end

%OLS loglik
olsloglik=-n*0.5*log(ssemat(1,1));
%Restricted logliks
restlogliks=[closeloglik;olsloglik];
%This forms the deviances
likratios=2*(maxlik-restlogliks);
%This gives the probability of observing a higher likelihood ratio under
%repeated sampling
prhigher=1-gammainc(likratios/2,1/2);
%augment bmax with alphamax
bmax=[bmax;alphamax];
%The signed root deviance is the square root of the deviance given the sign
%of the parameter
srds=sqrt(likratios).*sign(bmax);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [closeloglik, closealpha]=closeopt(ssemat, c, n)

%extracts eo'eo,ed'eo,ed'ed
eoteo=ssemat(1,1);
edteo=ssemat(1,2);
edted=ssemat(2,2);

%defining intermediate variables
a=(2*c-n)*edted;
b=(n-4*c)*edteo;
f=2*c*eoteo+n*edted;
g=-n*edteo;

a1=b/a;
a2=f/a;
a3=g/a;

%setting up cubic equation
Q=(3*a2-a1^2)/9;
R=(9*a1*a2-27*a3-2*a1^3)/54;
D=Q^3+R^2;

S=(R+sqrt(D))^(1/3);
T=(R-sqrt(D))^(1/3);

%possible solutions to cubic equation
alpha0=S+T-(1/3)*a1;
alpha1=-(1/2)*(S+T)-(1/3)*a1+(1/2)*i*sqrt(3)*(S-T);
alpha2=-(1/2)*(S+T)-(1/3)*a1-(1/2)*i*sqrt(3)*(S-T);

%vector of three possible solutions
alphavec=[alpha0;alpha1;alpha2];

%the optimal solution will fall between -1,1
alphaind=logical((alphavec>-1).*(alphavec<1));

%extracting optimal solution
closealpha=alphavec(alphaind);

%maximum liklihood value (up to a constant)
closeloglik=c*log(1-closealpha^2)-(n/2)*log(eoteo-2*closealpha*edteo+(closealpha^2)*edted);

