%This example program uses the data from:
%Pace, R. Kelley, and Ronald Barry, Quick Computation of Regressions with a Spatially Autoregressive Dependent Variable,
%Geographical Analysis, Volume 29, Number 3, July 1997, p. 232-247. 

%This program computes nearest neighbor spatial weight matrices based upon the 1st and 2nd order Delaunay triangle neighbors
%If one requests more neighbors than found in the first and second order neighbors for a particular observation, the neighbors beyond these will have
%zero elements in the individual neighbor weight matrices. Hence, this is a set of nearest neighbors restricted to the 1st and 2nd order Delaunay neighbors.

%Important parameters include the number of neighbors, m, and rho, the decline of the weight of the jth neighbor with j

%Written by Kelley Pace, www.spatial-statistics.com, 6/23/97 and revised
%12/26/02

%some standard clearing, and setting display format
clear all;
clear global;
close all;
format bank;

%loading election data
load y;
load xsub;
load xcoord;
load ycoord;

%number of observations
n=length(ycoord);

%Intercept variable
o=ones(n,1);

%standardizing dimensions
avxsub=mean(xsub);
stdxsub=std(xsub);
xsub=(xsub-avxsub(o,:))./stdxsub(o,:);

%standardizing space -- need to use common scaling
space_coord_scale=max(std(xcoord),std(ycoord));
space_coords=[(xcoord-mean(xcoord))/space_coord_scale (ycoord-mean(ycoord))/space_coord_scale ];


%nnmax is the maximum number of neighbors desired
nnmax=30;
%can be from 1 to nnmax
rho=0.95;
%m, the number of neighbors used in the actual weight matrix, can be from 1 to nnmax
m=30;

%some possible weights given the xsub dimensions
wvec=[  0.02 -0.02 0]';
%number of wvec parameteers
witer=length(wvec);
%a unit vector to match wvec
ovec=ones(witer,1);

%forms all possible weights (cases) given to the dimensions
%2 parms
%desmat=[kron(wvec,ovec) kron(ovec,wvec)];
%3 parms
%desmat=[kron(kron(wvec,ovec),ovec) kron(kron(ovec,wvec), ovec) kron(kron(ovec,ovec),wvec)];
%4 parms
desmat=[kron(kron(kron(wvec,ovec),ovec),ovec) kron(kron(kron(ovec,wvec), ovec),ovec) kron(kron(kron(ovec,ovec),wvec),ovec) kron(kron(kron(ovec,ovec),ovec),wvec) ];%3 parms


%Note, the pure spatial case is at the end for comparison (last row of
%desmat is all zeros.

%number of cases
[diter,ndparms]=size(desmat);

%looping over the cases
tic
%for i=1:diter
for i=1:1
    
%vector of dimension weighting parameters
desvec=desmat(i,:);
%forming matrix to weight dimensions (columns) in xsub
dwmat=diag(desvec);

%spatial coordinates, weighted non-spatial coordinates
coords=[space_coords xsub*dwmat];

%Computes nearest neighbors in all dimensions supplied
fneighbors_multi2(coords, nnmax);

%Symmetricizes and weights nearest neighbors
[wswnn,wwsnn,wmatnn]=fsym_neighbors2(m, rho);
%wswnn is a symmetric weight matrix
%wwsnn is a row-stochastic assymetric weight matrix
%wmatnn is the diagonal weighting matrix used in making the above

%forms x with spatial lags
x=[xsub wswnn*xsub o];

%Estimating MESS-AR
[bmax, srds, prhigher, emax, maxlik]=fmess_ar2(x, y, wswnn);

%recording weight parameters, spatial parameter, likelihood
rez(i,:)=[desmat(i,:) bmax(end) maxlik];

rez(i,:)

end;
toc

%maximum profile likelihood, and which row it is in
[maxlik,maxindc]=max(rez(:,end));

%results for multidimensional
extended_space=rez(maxindc,:);
%pure space ones -- since last iteration was pure spatial (all 0s in last
%row of desmat)
pure_space=rez(end,:);

%putting in optimum weights
%vector of dimension weighting parameters
desvecopt=extended_space(1:ndparms);
%forming matrix to weight dimensions (columns) in xsub
dwmatopt=diag(desvecopt);
coords=[space_coords xsub*dwmatopt];

%Computes nearest neighbors in all dimensions supplied
fneighbors_multi2(coords, nnmax)

%Symmetricizes and weights NNs 
[wswnn,wwsnn,wmatnn]=fsym_neighbors2(m, rho);

%forms x with spatial lags
x=[xsub wswnn*xsub o];

%Estimating MESS-AR for the multdimensional case
[bmaxm, srdsm, prhigherm, emaxm, maxlikm]=fmess_ar2(x, y, wswnn);


%printing using James LeSage's mprint, www.spatial-econometrics.com

disp(' ')
disp('Testing for the signficance of additional dimensions giving possible parameter values of -0.02, 0.00, and 0.02 to each extended dimension.')
disp('If the likelihood is unimodal in the extended dimensions, then it should result in an estimate other than 0, since both small positive')
disp('and negative values are possible. In other words, estimates at the small parameter values should inform about existence and ')
disp('direction of the multidimensional dependence.')
disp(' ')

info.width=132;
info.rnames=strvcat('Variables ' , 'Voting Pop ','Education ','Home Ownership ','Income ','Lag Voting Pop','Lag Education','Lag Home Ownership  ','Lag Income', 'Intercept','Alpha');
info.cnames=strvcat('b 2-D ', 'SRD 2-D ', 'b 6-D ', 'SRD 6-D ');

%pure spatial was last case -- can just use values from the loop
mprint([bmax srds bmaxm srdsm],info)

disp(' ')

info0.rnames=strvcat(' ','parameters');
info0.cnames=strvcat('Number of neighbors ',' Geometric parameter');
mprint([m rho],info0)
disp(' ')

info1.rnames=strvcat('Spaces','2-D ', '6-D ');
info1.cnames=strvcat('Pop-D','Edu-D', 'Home-D','Income-D  ',' Mess AR parameter ','Loglik');
mprint([pure_space;extended_space],info1);

disp(' ')
disp('The multidimensional clearly has a better likelihood, and materially affects some parameter estimates and SRDs.')
