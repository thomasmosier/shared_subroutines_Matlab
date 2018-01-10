function [doubly_d]=fdoubly2(d)
%
%[doubly_d]=fdoubly2(d)
%
%This function takes a non-negative, symmetric weight matrix and makes it doubly stochastic
%(both rows and columns sum to 1).
%
%INPUT:
%
%d is a n by n, non-negative, and symmetric weight matrix.
%
%OUTPUT:
%
%doubly_d is a n by n symmetric and doubly stochastic weight matrix.
%
%NOTES:
%
%If a matrix is overly sparse, this function may not converge. Also, this function can be somewhat slow for large matrices.
%You can adjust the maximum iterations or the convergence criterion in cond1
%to improve accuracy at the cost of performance and vice-versa. The
%function may work for some asymmetric matrices, but this should not be
%expected.
%
%If you use this function, please cite:
%
%%Pace, R. Kelley, and James P. LeSage, 
%"Semiparametric Maximum Likelihood Estimates of Spatial Dependence," 
%Geographical Analysis, Volume 34, Number 1, January 2002, p. 75-90.
%
%Kelley Pace, www.spatial-statistics.com, 12/25/02

%Starts by finding one of the dimensions of the n by n matrix d
n=size(d,1);

%initializes loop
i=0;
allcond=1;
ww=ones(n,1);
%this block of code alternates between normalizing the weight matrix for
%the row sums and the column sums. For symmetric matrices one could use
%either one instead.
while (allcond);
   i=i+1;
   wr=(1./sqrt(sum(d')))';%row sum normalization
   wc=(1./sqrt(sum(d)))';%column sum normalization
   
   if mod(i,2)>0;%uses wr or wc depending upon even or odd i
  w=wr;
else;
  w=wc; 
 end

%forms a diagonal matrix
wmat=spdiags(w,0,n,n);
 %normalizes by one of the sums
d=wmat*d*wmat;

%tests for convergence, or stops after 1000 loops
cond1=(max(abs(wr-1))+max(abs(wc-1)))>0.01;
cond2=(i<1000);
allcond=cond1*cond2;

end;

%prints messaage if convergence isn't reached before hitting maximum
%iterations
if i>999;
   ['Warning - May not have converged - iteration number, max,min of wr,wc below']
   ['These should not be off by too much (e.g., 0.1)']
   [i max(wr) min(wr) max(wc) min(wc)]
   end;
   
%returns answer
doubly_d=d;
