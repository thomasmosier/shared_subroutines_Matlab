%demonstrates sequence of log-determinant ideas using Chebyshev approximation.
%
%Written by Kelley Pace, www.spatial-statistics.com, 1/01/03.

clear all
close all
clear global
format bank
clc

n=5000;
xcoord=rand(n,1);
ycoord=rand(n,1);
[wswdel,wwsdel,wmatdel]=fdelw2(xcoord,ycoord);
d=wswdel;


tic;
 [lowerbounds, chebyshevest, upperbounds, alphafine]=fdet_chebyshev_seq2(d);
 chebyshev_sequence_time=toc;
 
 
plot(alphafine, [lowerbounds(:,n) chebyshevest(:,n) upperbounds(:,n)])
 xlabel('\alpha')
 ylabel('Log-determinant')
 title('Log-determinant versus spatial dependence parameter for all n observations')
 legend('Taylor lower bound','Chebyshev approximation','Taylor upper bound',3)
 
 figure
 
 nsub=round(n/2);
plot(alphafine, [lowerbounds(:,nsub) chebyshevest(:,nsub) upperbounds(:,nsub)])
 xlabel('\alpha')
 ylabel('Log-determinant')
 title('Log-determinant versus spatial dependence parameter for n/2 observations')
 legend('Taylor lower bound','Chebyshev approximation','Taylor upper bound',3)
 
 disp('As expected, the log-determinant plot for all observations has a similar appearance to')
 disp('the log-determinant plot for half the observations (but the magnitudes are roughly different')
 disp('by a factor of 2 as well')
 disp(' ')
 
  
p=symrcm(d);%symmetric reverse Cuthill-McKee is a local ordering
dlocal=d(p,p);%reordering to make local
figure
spy(dlocal)
xlabel('Observation Number')
ylabel('Observation Number')
title('Local Permutation Non-zero Weight Pattern for Delaunay')
disp('One can see the local nature of the spatial dependence in the reordered weight matrix')
disp(' ')


[lowerbounds, chebyshevest, upperbounds, alphafine]=fdet_chebyshev_seq2(dlocal);
 
 na=length(alphafine);%largest dependence parameter 
 nasub=round(na/2);%medium dependence parameter
 
 figure
 plot(chebyshevest([nasub na],:)')
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
disp('One can see the random nature of the spatial dependence in the reordered weight matrix')
disp(' ')

 [lowerbounds, chebyshevest, upperbounds, alphafine]=fdet_chebyshev_seq2(drand);
 
 figure
 plot(chebyshevest([nasub na],:)')
 legend(['alpha=' num2str(alphafine(nasub))],['alpha=' num2str(alphafine(na))],3)
 title('Log-determinant versus subsample size for two spatial dependence parameters (random ordering)')
 ylabel('Log-determinant')
 xlabel('sub-sample size')
 
 disp('In contrast, the random ordering shows a logarithmic pattern with respect to adding more observations,')
 disp('especially for large spatial parameters. This case corresponds to random sampling or infill asymptotics.')
 disp('Essentially, this increases the density of observations per area. Intuitively, there is less information')
 disp('being added for each additional observation and the logarithmic pattern shows the correct penalty.')
 disp('If the degree of spatial dependence is high, this should increase the logarithmic pattern as shown in')
 disp('in the plot with the larger dependence parameter.')
 disp(' ')
 disp('It takes little time to compute these sequences using the Chebyshev log-determinant approximation.')
 chebyshev_sequence_time
 disp('seconds')