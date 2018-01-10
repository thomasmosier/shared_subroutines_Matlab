function [Coefficients, S_err, XTXI, R_sq, F_val, Coef_stats, Y_hat, residuals, covariance] ...
          = mregress(Y, X, INTCPT)

      
%%COMMENTS FROM FILE EXCHANGE:
% ===========================================
% % Solve for the regression coefficients using ordinary least-squares 
% % ------------------------------------------------------------------
%    if (INTCPT == 1) 
%      X = [ ones(n,1) X ] ; 
%    end 
%   
%     XTXI = inv(X' * X); 
%     Coefficients = XTXI * X' * Y ;
% ===========================================
% See that the normal equations are used, as well as (shudder) inv.
% This is the WRONG way to solve this problem. It is a mistake that many 
%novices to regression make. It is a mistake that many people who have done 
%regression for years make, simply because they have never been told they 
%are doing it the wrong way. It is a mistake that many authors of 
%statistics texts make, simply because nobody has ever told them either.
% The normal equations, are
%    C = inv(X' * X) * X' * Y
% are ONE way to solve the over-determined problem
%    X*C = Y
% They are a POOR solution. First of all, this formula uses inv. Somewhat 
%better would be to write it as
%    C = (X' * X) \ (X' * Y)
% The use of backslash here avoids the use of inv. But it still fails 
%because the normal equations are still involved, just in another form.
% Look at a better solution.
%    C = X \ Y;
% This uses the built-in backslash to solve the problem. How would 
%backslash solve it internally? (I'll assume that no pivoting is done to 
%make things simple.) Suppose we choose to form a QR factorization of X.
%    X = Q*R
% Here, Q is an orthogonal matrix, R is upper triangular. (Again, I'll 
%ignore any pivoting considerations.) Our problem is now to solve
%   Q*R*C = Y
% Since Q is orthogonal, then we can left multiply by the transpose of Q, 
%to reduce the problem to
%   R*C = Q' * Y
% Remember, Q was orthogonal, so that (Q' * Q) was an identity matrix.
% We now solve for C as
%   C = R \ (Q' * Y)
% Backslash is smart here. It knows that R is upper triangular, so 
%backslash will do the solution as a back-substitution, an O(n^2) operation. 
%No explicit inverse is required at all.
% Again, backslash does all of this internally. But it is a far better 
%solution than the normal equations. What happens in the normal equations 
%to make that solution poor? The problem is the condition number of the 
%matrix. If your matrix is at all poorly conditioned (something that is 
%very common in multiple regression problems) then when you form (X'*X), 
%the condition number of the matrix will be SQUARED. So the normal 
%equations are a terribly poor way to solve the system posed.
% We can get the condition number of a matrix from the cond function. It is 
%defined as the ratio of the largest to the smallest singular value of the 
%matrix. When this number is large, then you should expect to loose 
%accuracy in the solution of your equations. Look at a simple, random 
%matrix. (By the way, the condition number will be small here. But most 
%regression problems have a much larger condition number.)
% X = rand(10,5); 
% svd(X) 
% ans = 
%           3.70360600537412 
%           1.16305761722647 
%          0.819599668934314 
%           0.77149157109882 
%          0.359840732709897
% See that the largest and smallest singular values would have a ratio of a wee bit over 10.
% cond(X) 
% ans = 
%           10.2923478881418
% What happens when I form X'*X?
% cond(X'*X) 
% ans = 
%           105.932425050538
% Yes. As expected, the condition number was squared. Try this again with a
%more typical regression problem. Here, just a simple 10th order polynomial in one variable.
% x = rand(20,1); 
% X = bsxfun(@power,x,[10 9 8 7 6 5 4 3 2 1 0]);
% cond(X) 
% ans = 
%           32218698.6067888
% cond(X'*X) 
% ans = 
%       1.03682639151296e+15
% See that the condition number of (X' * X) is roughly 1e15. In fact, if we 
%use the normal equations to solve this problem, the coefficients will have 
%few correct digits in them, because (X' * X) is nearly a numerically singular matrix.
% Now try an 11th order polynomial. Make some random data too. Since the 
%actual numbers don't matter, just see if there is a difference in the solution.
% X = bsxfun(@power,x,[11 10 9 8 7 6 5 4 3 2 1 0]); 
% Y = rand(20,1);
% X\Y 
% ans = 
%           213864.116016404 
%          -1151194.24570325 
%           2686678.25926054 
%          -3564335.59383269 
%           2961378.38376311 
%          -1600304.51887482 
%           565978.359556758 
%           -128793.21499217 
%           18153.6837583098 
%          -1482.77560777251 
%           58.7870441102109 
%         0.0478989431356096
% Now compare that to the normal equations. See that inv complains here, because the system is now effectively singular.
% inv(X'*X)*X'*Y 
% Warning: Matrix is close to singular or badly scaled. 
%          Results may be inaccurate. RCOND = 8.102118e-18. 
% ans = 
%           223829.921190437 
%          -1205573.39200913 
%           2815035.41413421 
%          -3736079.41837729 
%           3104808.97033064 
%          -1677906.49043969 
%           593313.714534699 
%          -134931.578398688 
%           18988.6654796978 
%          -1545.18491246855 
%            60.905800953616 
%         0.0268631240254516
% See that many of the resulting coefficients do not even agree in a single digit of the result. The problem is when we formed (X' * X).


% MREGRESS  Performs multiple linear regression analysis of X (independent) on Y (dependent).
%
%           Usage:
%
%           [Coefficients, S_err, XTXI, R_sq, F_val, Coef_stats, Y_hat, residuals, covariance] ... 
%                                                     = mregress(Y, X, INTCPT)
%
%	    INTCPT = 1; include a y-intercept in the model
%           INTCPT = 0; DO NOT include a y-intercept in the model
%
%           Returns:
%
%           Coefficients         - Regression coefficients
%           S_err                - Standard error of estimate
%           XTXI                 - inverse of X' * X
%           R_sq                 - R-squared
%           F_val                - F-value for the regression and siginificance level (p-value for F)
%           Coef_stats           - Coefficients with their standard deviations, T-values, and p-values
%           Y_hat                - Fitted values
%           residuals            - Residuals
%	    covariance           - Covariance matrix ( XTXI * S_err^2 )
%
%
% G. Anthony Reina
% Motor Control Lab
% The Neurosciences Institute
% Created: 4 Aug 1998
% Last Update: 10/8/1998 by GAR
%
% Please note that for the case when the intercept of the model equals zero, the
% definition of R-squared and the F-statistic change mathematically. For a linear
% model containing the y-intercept, R-squared refers to the amount of variance around
% the mean of the dependent variable (y) which is explained by the variance around the
% mean of the independent variables (x). For a linear model NOT containing the
% y-intercept, R-squared measures the amount of variance around ZERO of the dependent
% variable which is explained by the variance around ZERO of the independent variable.
% If the same equation for R-squared is used for both with and without a y-intercept
% (namely R-squared = [Sum of Squares of the Regression] / [ Total sum of the squares]),
% then R-squared may be a NEGATIVE value for some data. For this reason, 
% this subroutine will calculate R-squared using the total un-corrected sum of the
% squares. In effect, this approach avoids negative R-squares but may lack any
% meaningful interpretation for the "goodness-of-fit" in the model. It has been
% suggested by some texts that a more useful approach is to always use the case
% where y-intercept is included in the model. However, as with all statistical
% analyses, it is prudent to simply be aware of what your data "looks" like and
% what your statistical tools are actually measuring in order to generate a useful
% analysis.
%
% For further reading on regression through the origin (i.e. without a y-intercept),
% please refer to:
%
%    Neter J, Kutner MH, Nachtsheim CJ, and Wasserman W. "Applied Linear 
%             Statistical Models" 4th ed. Irwin publishing (Boston, 1996), pp 159-163.
%
%    Myers R, "Classical and Modern Regression with Applications" Duxbury Press
%             (Boston, 1986), p. 30.
%
%

if (nargin < 2)
   error('mregress requires at least 2 input variables. Type ''help mregress''.');
end

if (nargin == 2) 
   INTCPT = 0;
end

% Check that independent (X) and dependent (Y) data have compatible dimensions
% ----------------------------------------------------------------------------
[n_x, k] = size(X);
[n_y,columns] = size(Y);

if n_x ~= n_y, 
    error('The number of rows in Y must equal the number of rows in X.'); 
end 

if columns ~= 1, 
    error('Y must be a vector, not a matrix'); 
end

n = n_x;

%  Solve for the regression coefficients using ordinary least-squares
%  ------------------------------------------------------------------

   if (INTCPT == 1)
     X = [ ones(n,1) X ]   ; 
   end
 
    XTXI = inv(X' * X);
    
    Coefficients = XTXI * X' * Y ;
    
    Coefficients = X\Y; 

%  Calculate the fitted regression values
%  --------------------------------------

    Y_hat = X * Coefficients;

%  Calculate R-squared
%  -------------------
% The calculation used for R-squared and the F-statistic here are based
% on the total, un-corrected sum of the squares as describe by Neter and
% Myers. Note that the meaning of R-squared changes for the case of
% regression without a y-intercept. This approach will yield the same
% results as SysStat, SAS, SPSS and BMDP but will differ from that of
% Excel, Quattro Pro, and the MATLAB regress.m function (for the case of
% no y-intercept in the model -- all packages including this one will
% agree for the case of linear regression with a
% y-intercept). Essentially, it is wise to find a way to
% keep the y-intercept (even if it is near zero) in the model to analyze
% it in a meaningful way that everyone can understand.

if (INTCPT == 1)

   RSS = norm(Y_hat - mean(Y))^2;   % Regression sum of squares.
   TSS = norm(Y - mean(Y))^2;       % Total sum of squares (regression plus residual).
   R_sq = RSS / TSS;                % R-squared statistic.

else

   RSS = norm(Y_hat)^2;             % Regression sum of squares.
   TSS = norm(Y)^2;                 % Total, un-corrected sum of squares.
   R_sq = RSS / TSS;                % R-squared statistic.

end

% $$$ % Alternative calculation of R-squared
% $$$ % ====================================
% $$$ % The follwing equation is from Judge G, et al. "An Introduction to the theory
% $$$ % and practice of econometrics", New York : Wiley, 1982. It is the
% $$$ % squared (Pearson) correlation coefficient between the predicted and
% $$$ % dependent variables. It is the same equation regardless of whether an
% $$$ % intercept is included in the model; however, it may yield a negative
% $$$ % R-squared for a particularily bad fit.
% $$$ covariance_Y_hat_and_Y = (Y_hat - mean(Y_hat))' * (Y - mean(Y));
% $$$ covariance_Y_hat_and_Y_hat = (Y_hat - mean(Y_hat))' * (Y_hat - mean(Y_hat));
% $$$ covariance_Y_and_Y = (Y - mean(Y))' * (Y - mean(Y));
% $$$ R_sq = (covariance_Y_hat_and_Y / covariance_Y_hat_and_Y_hat) * ...
% $$$        (covariance_Y_hat_and_Y / covariance_Y_and_Y);


%  Calculate residuals and standard error
%  --------------------------------------

    residuals = Y - Y_hat;

    if (INTCPT == 1)
       S_err = sqrt(residuals' * residuals / (n - k - 1) );
    else
       S_err = sqrt(residuals' * residuals / (n - k) );
    end

%  Calculate the standard deviation and t-values for the regression coefficients
%  -----------------------------------------------------------------------------

    covariance = XTXI .* S_err^2;
    
    C = sqrt(diag(covariance, 0));

    % (n.b. Need to perform a 2-tailed t-test)
    % ****************************************
    if (INTCPT == 1)
      p_value = 2 * (1 - tcdf(abs(Coefficients./C), (n - (k + 1))));
    else
      p_value = 2 * (1 - tcdf(abs(Coefficients./C), (n - k)));
    end

    Coef_stats = [ Coefficients, C, (Coefficients./C), p_value];


% Estimator of error variance.
% ----------------------------

if (INTCPT == 1) 

     SSR_residuals = norm(Y - Y_hat)^2;
     TSS = norm(Y - mean(Y))^2;     % Total sum of squares (regression plus residual).

     F_val = (TSS - SSR_residuals) / k / ( SSR_residuals / (n - (k + 1)));

     F_val = [F_val (1 - fcdf(F_val, k, (n - (k + 1)))) ];

else

     SSR_residuals = norm(Y - Y_hat)^2;
     TSS = norm(Y)^2;                % Total sum of squares (regression plus residual).

     F_val = (TSS - SSR_residuals) / k / ( SSR_residuals / (n-k));

     F_val = [F_val (1 - fcdf(F_val, k, (n - k))) ];

end


