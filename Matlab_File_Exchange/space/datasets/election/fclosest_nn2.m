function [nnlist]=fclosest_nn2(xcoord, ycoord)
%
%[nnlist]=fclosest_nn2(xcoord, ycoord)
%
%Finds the closest neighbor to each observation
%
%INPUT:
%
%The two n by 1 vectors xcoord, ycoord which represent locational
%coordinates.
%
%OUTPUT:
%
%The n by 1 vector nnlist containing indices to the nearest neighbor
%
%NOTES:
%
%Closest neighbor dependence was studied in:
%
%Pace, R. Kelley, and Dongya Zou, 
%“Closed-Form Maximum Likelihood Estimates of Nearest Neighbor Spatial Dependence,” 
%Geographical Analysis, Volume 32, Number 2, April 2000, p. 154-172.
%
%Written by Kelley Pace, www.spatial-statistics.com, 8/98, revised
%12/26/02.

%computes Delaunay triangles - returns a vector with 3 columns. Each column contains
%a reference to the row of the n by 2 vector of xcoord,ycoord. Hence, the three columns
%define a triangle and hence 6 possible pairs of connected points
tri=delaunay(xcoord, ycoord);

%this stacks the six possible pairs of connected points
bigm=[tri(:,[1 2]); tri(:,[1 3]); tri(:,[2 3]); tri(:,[2 1]); tri(:,[3 1]); tri(:,[3 2])] ;

%tri is not needed anymore, so to conserve memory:
clear tri;

%this computes squared Euclidean distance between the points in each pair.
%Ordinally, the squared distance is the same as the distance itself.
d=(xcoord(bigm(:,1))-xcoord(bigm(:,2))).^2+(ycoord(bigm(:,1))-ycoord(bigm(:,2))).^2;

%this sorts by row and by distance
pair_plus_d=sortrows([bigm d],[1 3]);

%bigm_d is no longer needed, so to conserve memory:
clear bigm d;

%this finds the beginning of each block of possible neighbors to observation i
difind=[1;diff(pair_plus_d(:,1))];

%this extacts the beginning row of each block (which has the smallest distance
%after the sort
nnlist=pair_plus_d(logical(difind),2);





