function [bmax, srds, prhigher, emax, maxlik]=fols2(x, y)
%
%[ bmax, srds, prhigher, emax, maxlik]=fols2(x, y)
%
%This OLS routine mimics the format of the spatial routines, and thus computes
%signed root deviance instead of t-ratios.
%
%INPUT:
%
%x represents n observations on k independent variables.
%
%y represents n observations on the dependent variable.
%
%OUTPUT:
%
%logliks is a (k+1) by 1 vector of optimum likelihoods, where the first
%element contains the optimum likelihood for the global model and the
%subsequent elements contain the optimum likelihoods for the delete-1
%submodels.
%
%emax is the n by 1 vector of residuals
%
%bmax is the k+1 by 1 vector of parameter estimates with the last element being a 0. 
%This aids with other routines
%
%SRDS is the k+1 by 1 vector of signed root deviances with the last element being a 0.
%
%prhigher is the k+1 by 1 vector of obtaining a higher SRD in repeated
%sampling. The last element is a 0.
%
%NOTES:
%
%The vectors bmax, SRDS, prhigher have been padded with a trailing 0 to make these match the spatial routines.
%
%Some of the tricks to compute the likelihood ratios appear in:
%
%Pace, R. Kelley, and Ronald Barry, Quick Computation of Regressions with a Spatially Autoregressive Dependent Variable,
%Geographical Analysis, Volume 29, Number 3, July 1997, p. 232-247. 
%
%Written by Kelley Pace, 6/23/97, and revised 12/25/02.

[n,k]=size(x);

bmax=x\y;
emax=y-x*bmax;
sse=emax'*emax;
maxlikols=-(n/2)*log(sse);

%this begins the computations of the delete 1 variable subsets for the restricted log-likelihoods
%this allows one to compare the likelihood ratios under OLS with those of other procedures
xtx=x'*x;
cholxtx=chol(xtx);
eyer=eye(k);

bz=(cholxtx')\eyer;
bsqterm=((bmax.^2))./(sum(bz.^2)');

rloglik=-(n/2)*log(sse+bsqterm);
likratios=2*abs(maxlikols-rloglik);
prhigher=1-gammainc(likratios/2,1/2);

maxlik=[maxlikols; rloglik]';

%The signed root deviance is the square root of the deviance given the sign
%of the parameter. It has an interpretation similar to t-ratios.
srds=sqrt(likratios).*sign(bmax);
bmax=[bmax;0];
srds=[srds;0];
prhigher=[prhigher;0];

