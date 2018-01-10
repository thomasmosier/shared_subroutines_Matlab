function r2 = pearson_cf(cf,x,y)

%Designed to use coeffients output from polyfit
yHat = polyval(cf,x);


%Calculation of R-squared
% ====================================
% The follwing equation is from Judge G, et al. "An Introduction to the theory
% and practice of econometrics", New York : Wiley, 1982. It is the
% squared (Pearson) correlation coefficient between the predicted and
% dependent variables. It is the same equation regardless of whether an
% intercept is included in the model; however, it may yield a negative
% R-squared for a particularily bad fit.
covariance_yHat_and_Y = (yHat - mean(yHat))' * (y - mean(y));
covariance_yHat_and_yHat = (yHat - mean(yHat))' * (yHat - mean(yHat));
covariance_Y_and_Y = (y - mean(y))' * (y - mean(y));
r2 = (covariance_yHat_and_Y / covariance_yHat_and_yHat) * ...
    (covariance_yHat_and_Y / covariance_Y_and_Y);



% RSS = norm(yHat)^2;     % Regression sum of squares.
% TSS = norm(y)^2;        % Total, un-corrected sum of squares.
% r2 = RSS / TSS;      	% R-squared statistic.