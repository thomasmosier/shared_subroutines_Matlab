function output = quad_fit_grids(input, xIn, xOut)

if numel(xIn) ~= 3
    error('linFitGrid:xInnumel',['The input x vector has ' num2str(numel(xIn)) ' entries, but length must be 3.'])
end
szIn = size(input);
if szIn(1) ~= 3
  error('linFitGrid:yInNumel',['The input y array has ' num2str(szIn(1)) ' entries, but length of first indice must be 3.'])
end

%Use Lagrange polynomial:
%http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html
output = nan([numel(xOut), numel(input(1,:,1)), numel(input(1,1,:))], 'single');
for ii = 1 : numel(xOut)
    output(ii,:,:) = ...
          (((xOut(ii)-xIn(2))*(xOut(ii)-xIn(3)))/((xIn(1)-xIn(2))*(xIn(1)-xIn(3))))*squeeze(input(1,:,:)) ...
        + (((xOut(ii)-xIn(1))*(xOut(ii)-xIn(3)))/((xIn(2)-xIn(1))*(xIn(2)-xIn(3))))*squeeze(input(2,:,:)) ...
        + (((xOut(ii)-xIn(1))*(xOut(ii)-xIn(2)))/((xIn(3)-xIn(1))*(xIn(3)-xIn(2))))*squeeze(input(3,:,:));
end

output = squeeze(output);

% %USE EXACT-FIT QUADRATIC TO INTERPOLATE DAYS:
% %These fitting parameters are matrices the same size as the
% %climate grids
% aFit = ( (input(2,:,:) - input(1,:,:))*(xIn(1) - xIn(3)) + (input(3,:,:) - input(1,:,:))*(xIn(2) - xIn(1)) ) ...
%     / ( (xIn(1) - xIn(3))*(xIn(2)^2 - xIn(1)^2) + (xIn(2) - xIn(1))*(xIn(3)^2 - xIn(1)^2) );
% 
% bFit = ( (input(2,:,:) - input(1,:,:)) - aFit*(xIn(2)^2 - xIn(1)^2) ) / (xIn(2) - xIn(1));
% 
% cFit = (input(1,:,:) - aFit*xIn(1)^2 - bFit*xIn(1));
% 
% %Apply coeffifient matrices to desired output points:
% output = nan([numel(xOut), numel(input(1,:,1)), numel(input(1,1,:))], 'single');
% for ii = 1 : numel(xOut)
%     output(xOut, :,:) = aFit*xOut^2 + bFit*xOut + cFit;
% end