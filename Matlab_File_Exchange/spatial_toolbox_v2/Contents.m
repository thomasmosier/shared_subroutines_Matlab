% Spatial Statistics Toolbox 2.0
% Kelley Pace, www.spatial-statistics.com, 1/15/03 
% 
% Two dimensional weight matrices
% 
% fdelw2 - spatial contiguity weight matrix
% fneighbors2 - nearest neighbors weight matrix in two dimensional space
% fclosest_nn2 - closest neighbor weight matrix
% 
% 
% Multidimensional weight matrices
% 
% fneighbors_multi2 - nearest neighbors weight matrix in 2+ dimensional space
% 
% 
% Spatiotemporal weight matrices
% 
% fnn_ml2 - spatiotemporal nearest neighbors in matlab
% fnn_mex2 - spatiotemporal nearest neighbors in Fortran 90 Mex file
% fts - spatiotemporal weight matrix with geometrically declining weights
% movav2 - temporal weight matrix
% 
% 
% Weight matrix reweighting functions
% 
% fasym_neighbors2 - row-stochastic weight matrix with geometrically declining weights
% fsym_neighbors2 - symmetric weight matrix with geometrically declining weights
% 
%
% Weight matrix scaling functions
%
% fdoubly2 - converts symmetric weight matrix into doubly stochastic weight matrix
%
% 
% Log-determinant functions
% 
% fdet_interp2 - exact log-determinant computations
% fdet_mc2 - Monte Carlo log-determinant estimation with confidence limits
% fdet_chebyshev2 - Chebyshev log-determinant estimation with Taylor series bounds
% fdet_interp_seq2 - sequence of exact log-determinants
% fdet_chebyshev_seq2 - sequence of approximate log-determinants
% 
% 
% Statistical functions
% 
% fols2 - OLS with signed root deviance inference
% fclosest_ar2 - closest neighbor closed-form autoregressive estimation via maximum likelihood
% fmess_ar2 - matrix exponential closed-form autoregressive estimation via maximum likelihood
% fmix2 - mixed regressive spatially autoregressive model estimation via maximum likelihood
% fsar2 - SAR estimation via maximum likelihood
% fcar2 - CAR estimation via maximum likelihood
% fmess_car2 - CAR (or SAR) estimation via maximum likelihod of the matrix exponential spatial specification
% frecursive2 - recursive estimation (used in spatiotemporal functions)
% 
% 
% Local estimation
% 
% fwholesale2 - spatial autoregressive local estimation (SALE)
% 
% 
% Simulation functions
% 
% fmess_sim2 - simulates MESS random variables
% fsar_sim2 - simulates SAR random variables
% fcar_sim2 - simulates CAR random variables
% 
% 
% Utilities
% 
% fdotplot2 - plots a variable over space with symbol size increasing with the values of the variable
% fmin_eig2 - estimates minimum eigenvalue of a weight matrix
% fgrid_generate2 - generates optional grid of values of dependence parameter used in some routines
% latextab2 - prints table to latex table with decimal alignment and commas every three digits
% tabprint2 - creates structure for printing with decimal alignment and commas every three digits.
% tab2 - prints table to screen or file with title (optional) and many features. Derived from James LeSage's mprint




