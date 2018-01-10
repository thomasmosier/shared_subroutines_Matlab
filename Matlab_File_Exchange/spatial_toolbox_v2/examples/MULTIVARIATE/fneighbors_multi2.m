function fneighbors_multi2(coords, nnmax)
%
%function fneighbors_multi2(coords, nnmax)
%
%This function finds the nnmax nearest neighbors to each observation and stores it in a series of individual neighbor matrices.
%Users should set mmax to be equal or higher to the number of neighbors they will
%need in the future.  This routine depends upon the TSTOOL package, which
%must be installed for this routine to work. Go to
%http://www.physik3.gwdg.de/tstool/ to download the package.
%
%INPUT:
%
%The n by nd locational coordinate matrix coords, where nd represents the
%number of dimensions.
%
%The scalar nnmax which is the largest number of neighbors needed in future
%spatial regressions.
%
%OUTPUT:
%
%The disk file smats.mat which contains a sequence of nnmax sparse matrices each
%containing at most a single element 1 in each row. These are labeled s1,
%s2, s3, .... snnmax. The functions fsym_neighbors2 and fasym_neighbors2
%work with smat and its component files.
%
%NOTES:
%
%One can optimize the likelihood over neighbors and weighting of neighbors
%as in:
%
%Pace, R. Kelley and Ronald Barry, O.W. Gilley, C.F. Sirmans, 
%“A Method for Spatial-temporal Forecasting with an Application to Real Estate Prices,” 
%International Journal of Forecasting, Volume 16, Number 2, April-June 2000, p. 229-246.
%
%Written by Kelley Pace, www.spatial-statistics.com, 12/26/02. However, the
%real essence of this routine is a call to the TSTOOL commands. The TSTOOL package
%is found at  http://www.physik3.gwdg.de/tstool/ and is intended for non-linear time series, but has some nice routines.
%The TSTOOL functions are MEX files and must be installed for your platform. However, the TSTOOL package
%comes with the MEX functions precompiled for multiple platforms.
%
%At this time, the TSTOOL routine is faster and more general than the Delaunay-based routine
%fneighbors2. Note, it does require a very painless installation of the
%TSTOOL routines.
%

%This should alert users to the need for TSTOOL
test_tstool=( exist('nn_prepare')*exist('nn_search') )==0;

if test_tstool==1
    error('TSTOOL does not appear to be installed, or the path to it is not correct. Go to http://www.physik3.gwdg.de/tstool/ to get the package, and install TSTOOL.');
    return
end

%obtaining the number of observations and dimensions
[n,nd]=size(coords);

%preparing the data structure -- see TSTOOL documentation
atria=nn_prepare(coords);

%There are other options for these routines such as precsion which could
%affect execution times. See the TSTOOL documentation for details.

%finding the indices of the nearest neighbors and their distances -- from
%TSTOOL
[index, distance]=nn_search(coords, atria, coords, (nnmax+1));

%this next set of statements converts the indices into m individual weight
%matrices

rowseqs=index(:,1);%row sequence
vals1=ones(n,1);%value of 1 for non-zero elements
vals0=zeros(n,1);%some zero elements to force shape

for i=1:nnmax;

%starts with the second column and goes until the nnmax+1 column -- hence
%nnmax neighbors
colseqs=index(:, (i+1));

%places 1s at the neighbors in the adjacency matrix
z1=[rowseqs colseqs vals1];


%these next two statements insure the dimensions are right. By placing a
%diagonal of zero, it forces the matrix to be n by n. For this function,
%z2, and vals0 may be optional.

z2=[rowseqs rowseqs vals0];
z=[z1;z2];

%converts index into an individual sparse matrix
eval(['s' num2str(i) '=spconvert(z)' ';']);

end;

%saves the individual matrices into one variable, and places this on disk.
%This file can be used by other routines.
save smats s* ;


