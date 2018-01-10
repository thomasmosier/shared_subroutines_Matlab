function fdotplot2(xcoord, ycoord, v, quantiles, gweights)
%
%fdotplot2(xcoord, ycoord, v, quantiles, gweights)
%
%This function plots dots of different weight and color for the intervals
%defined by the quantile vector at the coordinates xcoord, ycoord.
%
%INPUT:
%
% xcoord, ycoord are n by 1 vectors representing the two-dimensional
% coordinates.
%
% v is a n by 1 vector. Usually v=f(xcoord, ycoord) where f(.) some
% function, possibly non-smooth.
%
% quantiles is a q element vector giving the quantiles (q>1). These range between
% 0 to 100. The q quantiles define q-1 intervals. There has to be at least
% 2 elements in the vector.
%
% gweights is a (q-1) element vector giving to weight given to dots from
% different quantile intervals.
%
%OUTPUT:
%
%A map of v versus the locational coordinates with differnent dot sizes and
%colors.
%
%NOTES:
%
%You can set different colors, and the repeatability of behavior by
%uncommenting  some of the commands below.
%
%
%Written by R. Kelley Pace, www.spatial-statistics.com, 12/22/2002.



% As a default, this will look at the intervals 0-25, 25-75, 75-100
if nargin<4
    quantiles=[0 25 75 100];
end;

% The number of intervals is one less the number of quantiles
n_intervals=length(quantiles)-1;

% As a default, this sets up dot weights starting at 8 points-squared and rising linearly
if nargin<5
    gweights=(8:8:(8*n_intervals));
end

% If you wish a different colormap, please change to taste.
cmat=colormap('prism');
%cmat=colormap('gray');

% If you want different behavior each invocation, uncomment:
%cmat=cmat(randperm(size(cmat,1)),:);

% If you wish black, uncomment below:
%cmat=zeros(100,3);

%number of points to plot
n=length(v);

% We sort the vector in increasing value and reorder the x,y coordinantes
% to match
[vs, vind]=sort(v);
xcoords=xcoord(vind);
ycoords=ycoord(vind);

% This gives a vector of the indices associated with a change in quantiles
nvec=round((n-1)*quantiles*0.01)+1;

% svec must have a positive minimum. I set it to 1.
svec=ones(n,1);

% I set the colors to white, in case someone misspecifies the quantiles.
% This way one can map values between any set of quantiles
c=ones(n,3);

%In the loop, I set the dot weights and colors by interval
for j=1:n_intervals;
    svec(nvec(j):nvec(j+1))=gweights(j);
    c(nvec(j):nvec(j+1),1)=cmat(j,1);
    c(nvec(j):nvec(j+1),2)=cmat(j,2);
    c(nvec(j):nvec(j+1),3)=cmat(j,3);
end

% I invoke the Matlab intrinsic command scatter. If you wish circles, erase
%  'filled'. You can change markers to something else as well.
scatter(xcoords,ycoords,svec,c,'filled')
