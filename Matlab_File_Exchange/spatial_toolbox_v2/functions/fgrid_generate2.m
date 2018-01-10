function [coarse_grid, fine_grid]=fgrid_generate2(ninterp, ntotal, abound)
%
%[coarse_grid, fine_grid]=fgrid_generate2(ninterp, ntotal, ebound)
%
%This function computes a fine and coarse grid of values for the spatial dependence parameter used in some
%of the routines.
%
%INPUT:
%
%The scalar coarse_grid gives the number of values of the spatial dependence parameter to use in computing the desired quantity.
%
%The scalar fine_grid gives the number of values of the spatial dependence parameter to use in finding the output answer. The
%difference among contiguous values of this grid govern the final precision of the answers.
%
%The scalar abound gives the lower bound of the dependence parameter. This has an effect only when the lower bound is
%negative. In the usual case, this would equal the reciprocal of the minimum eigenvalue of the spatial weight matrix.
%For all other supplied values of abound, the lower bound of the spatial dependence parameter is 0.
%
%OUTPUT:
%
%coarse_grid has ninterp+1 values, if the lower limit of the grid is 0, and 2*ninterp+1 values otherwise.
%
%fine_grid has ntotal+1 values, if the lower limit of the grid is 0, and 2*ntotal+1 values otherwise.
%
%NOTES:
%
%Written by Kelley Pace, www.spatial-statistics.com,12/31/02.

if ninterp<15;
    disp('ninterp may be too low for much accuracy')
end

ngsub=ninterp+1;
ntotals=ntotal+1;

maxeig=1;
maxalpha=maxeig-0.005;
gpos=linspace(0,maxalpha,ntotals)';
logdetestpos=log(1-0.9*gpos);
logpointspos=linspace(logdetestpos(1),logdetestpos(end),ngsub);
apointspos=spline(logdetestpos,gpos,logpointspos);
apointspos(1)=0;
coarse_grid=apointspos';
fine_grid=gpos;


if abound<0;
minalpha=abound+0.005;
mincomplement=-minalpha;
gneg=linspace(0,mincomplement,ntotals)';
logdetestneg=log(1+0.9*gneg/mincomplement);
logpointsneg=linspace(logdetestneg(1),logdetestneg(end),ngsub);
apointsneg=spline(logdetestneg,gneg,logpointsneg);
coarse_grid=union(-apointsneg(2:end),apointspos)';
fine_grid=union(-gneg,gpos);
end






