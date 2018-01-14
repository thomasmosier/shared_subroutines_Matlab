function output = Gaussian_filter(input, rad, varargin)

output = nan(size(input), 'single');

%The 'sigma' (kernel dissipation) can be specified using a variable input 
%argument. OTherwise it defaults to 1/3 of the radius.  
if~isempty(varargin(:))
    sigma = varargin{1};
else
    sigma = rad/3; 
end

%Make kernel:
xRad = (-rad : rad); 
yRad = (-rad : rad)';  %+/-  1 grid cell (half each side)
[xRadMesh, yRadMesh] = meshgrid(xRad,yRad);
f = exp(-xRadMesh.^2/(2*sigma^2) - yRadMesh.^2/(2*sigma^2));
f = f./sum(f(:));

dataMean = mean2d(input);

%Low pass filter:
%Use 'nanconv' so that nan values do not impact convolution
df = nanconv(input - dataMean, f, 'same');

%Filtered data:
output = df./conv2(ones(size(output), 'single'),f,'same') + dataMean;