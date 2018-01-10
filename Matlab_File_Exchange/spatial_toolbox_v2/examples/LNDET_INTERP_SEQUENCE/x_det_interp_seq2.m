%This script demos the sequence of log-determinant routine and shows how ordering of observations affects results.
%
%Kelley Pace, www.spatial-statistics.com, 1/01/03.

clear all
close all
clear global
format bank
clc

%loading county weight matrix
load d;

ndsub=1500;

d=d(1:ndsub,1:ndsub);%reducing the number of observations to speed running time for the example

tic;
[detvalz, alphafine]=fdet_interp_seq2(d, 1);
original_time=toc;
 
n=size(d,1);
 
plot(alphafine, detvalz(:,n) )
 xlabel('\alpha')
 ylabel('Log-determinant')
 title('Log-determinant versus spatial dependence parameter for all n observations')
 
figure
 
nsub=round(n/2);
plot(alphafine, detvalz(:,nsub) )
xlabel('\alpha')
ylabel('Log-determinant')
title('Log-determinant versus spatial dependence parameter for n/2 observations')
 
 disp('As expected, the log-determinant plot for all observations has a similar appearance to')
 disp('the log-determinant plot for half the observations (but the magnitudes are roughly different')
 disp('by a factor of 2 as well')
 disp(' ')
  
 figure
spy(d)
xlabel('Observation Number')
ylabel('Observation Number')
title('Original Non-zero Weight Pattern for Delaunay on County data')
disp('One can see somewhat local spatial dependence in the spy plot of the original weight matrix.')
disp(' ')
  
 na=length(alphafine);
 nasub=round(na/2);
 
  figure
 plot(detvalz([nasub na],:)')
 legend(['alpha=' num2str(alphafine(nasub))],['alpha=' num2str(alphafine(na))],3)
 title('Log-determinant versus subsample size for two spatial dependence parameters (original ordering)')
 ylabel('Log-determinant')
 xlabel('sub-sample size')
  
 disp('The orginal ordering, based upon grouping states, suggests somewhat linear behavior in the log-determinant')
 disp('as subsample size increases. However, this is bumpy due to the clustering of counties by state.')
 disp(' ')

p=symrcm(d);%symmetric reverse Cuthill-McKee is a local ordering
dlocal=d(p,p);%reordering to make local
figure
spy(dlocal)
xlabel('Observation Number')
ylabel('Observation Number')
title('Local Permutation Non-zero Weight Pattern for Delaunay')
disp('One can see the local nature of the spatial dependence in the spy plot of the reordered weight matrix')
disp(' ')

tic
[detvalz, alphafine]=fdet_interp_seq2(dlocal, 1);
local_time=toc;
 
figure
 plot(detvalz([nasub na],:)');
 legend(['alpha=' num2str(alphafine(nasub))],['alpha=' num2str(alphafine(na))],3)
 title('Log-determinant versus subsample size for two spatial dependence parameters (local ordering)')
 ylabel('Log-determinant')
 xlabel('sub-sample size')
 
 disp('The local ordering shows approximately linear behavior in the log-determinant, implying an almost')
 disp('constant cost or penalty as one adds an additional observation. This case corresponds to increasing')
 disp('domain asymptotics. That is, what happens as one increases the area of spatial data while collecting')
 disp('a constant density of observations.')
 disp(' ')
 
p=randperm(n);
drand=d(p,p);

figure
spy(drand)
xlabel('Observation Number')
ylabel('Observation Number')
title('Random Permutation Non-zero Weight Pattern for Delaunay')
disp('One can see the random nature of the spatial dependence in the spy plot of the randomly reordered matrix')
disp(' ')
tic;
 [detvalz, alphafine]=fdet_interp_seq2(drand, 1);
 random_time=toc;
 
 figure
 plot(detvalz([nasub na],:)')
 legend(['alpha=' num2str(alphafine(nasub))],['alpha=' num2str(alphafine(na))],3)
 title('Log-determinant versus subsample size for two spatial dependence parameters (random ordering)')
 ylabel('Log-determinant')
 xlabel('sub-sample size')
 
 disp('In contrast, to the local  ordering, the random ordering shows a logarithmic pattern, ')
 disp('especially for large spatial parameters. This case corresponds to random sampling or infill asymptotics.')
 disp('Essentially, this increases the density of observations per area. Intuitively, there is less information')
 disp('being added for each additional observation and the logarithmic pattern shows the correct penalty.')
 disp('If the degree of spatial dependence is high, this should increase the logarithmic pattern as shown in')
 disp('in the plot with the larger dependence parameter.')
 disp(' ')
 
 disp(' ')
  disp('Sequence of Log-determinant computing times as a function of orderings. ')
 info.cnames=strvcat('orignal ordering','local ordering','random ordering');
 info.rnames=strvcat('ORDERINGS: ','Time in seconds');
 %using James LeSage's mprint function (www.spatial-econometrics.com)
  tmat=[original_time local_time random_time];
  mprint(tmat,info)
 disp(' ')
 disp('It takes little time to compute these sequences for local orderings, but the time rises with the random ordering.')
 disp(' ')
 disp('If you ran this again with a larger value of ndsub, you would see a very steep increase in computing time')
 disp('of the log-determinant with random ordering.')