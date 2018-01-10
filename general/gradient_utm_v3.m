function [varargout] = gradient_utm_v3(data, yGrid, xGrid, utmZone, varargin)

%Extension of Matlab's gradient function for use with points described in
%UTM coordinates
%Special case where point straddles UTM zones is programmed for (uses 'modal' UTM zone for local point).


%Explanation from Matlab's gradient function:
%gradient calculates the central difference for interior data points. 
%For example, consider a matrix with unit-spaced data, A, that has 
%horizontal gradient G = gradient(A). The interior gradient values, G(:,j), are
% 
% G(:,j) = 0.5*(A(:,j+1) - A(:,j-1));
% The subscript j varies between 2 and N-1, with N = size(A,2).
% 
% gradient calculates values along the edges of the matrix with single-sided differences:
% 
% G(:,1) = A(:,2) - A(:,1);
% G(:,N) = A(:,N) - A(:,N-1);
% If you specify the point spacing, then gradient scales the differences 
%appropriately. If you specify two or more outputs, then the function also 
%calculates differences along other dimensions in a similar manner. 
%Unlike the diff function, gradient returns an array with the same number of elements as the input.


%See output arguments at bottom:
blPer = 0;
blMag = 0;
if ~isempty(varargin)
   for ii = 1 : numel(varargin(:)) 
       if strcmpi(varargin{ii}, 'per') %Calculate slope/change in quantity in fractional terms
           blPer = 1;
       elseif strcmpi(varargin{ii}, 'gradmag') %Output magnitude (swrt sum of squares)
           blMag = 1;
%        elseif regexpbl(varargin{ii}, 'delta') %Output magnitude (swrt sum of squares)
%            blDelta = 1;
       else
           error('gradientUtm:varargin', [varargin{ii} ' is an unknown option.']);
       end
   end
end

%Remove singleton dimensions
data = squeeze(data);

if ~ismatrix(data)
    error('gradient_utm:notMatrix','The input must be a matrix.');
end

%Initialize outputs:
xGrad = nan(size(data));
yGrad = xGrad;
deltaDataX = xGrad;
deltaDataY = xGrad;
deltaDistX = xGrad;
deltaDistY = xGrad;

%Calculate gradients at each point:
for ii = 1 : numel(data(:,1))
    for jj = 1 : numel(data(1,:))
        
        %Define latitude indices to use:
        if ii == 1
            indLat = [ii, ii + 1];
        elseif ii == numel(data(:,1))
            indLat = [ii-1, ii];
        else
            indLat = [ii-1, ii, ii + 1];
        end
        
        %Define longitude indices to use:
        if jj == 1
            indLon = [jj, jj + 1];
        elseif jj == numel(data(1,:))
            indLon = [jj-1, jj];
        else
            indLon = [jj-1, jj, jj + 1];
        end
        
        %Calculate x distance:
        if all(strcmpi(utmZone(ii, indLon), utmZone(ii, jj))) %Same UTM zone
            deltaDistXTemp = diff(xGrid(ii, indLon));
        else %Different UTM zones
            zoneUse = cell_mode(utmZone(ii, indLon));
            
            [xForced, ~] = utm_same_zone(xGrid(ii, indLon), yGrid(ii, indLon), utmZone(ii, indLon), zoneUse);
            deltaDistXTemp = diff(xForced);
        end
        
        %Calculate y distance:
        if all(strcmpi(utmZone(indLat, jj), utmZone(ii, jj))) %Same UTM zone
            deltaDistYTemp = diff(yGrid(indLat, jj));
        else %Different UTM zones
            zoneUse = cell_mode(utmZone(indLat, jj));
            
            [~, yForced] = utm_same_zone(xGrid(indLat, jj), yGrid(indLat, jj), utmZone(indLat, jj), zoneUse);
            deltaDistYTemp = diff(yForced);
        end
        
        deltaDataXTemp =  data(ii, indLon(end)) - data(ii,  indLon(1));
        deltaDataYTemp = (data(indLat(end), jj) - data(indLat(1), jj));
        deltaDistX(ii,jj) = mean(deltaDistXTemp);
        deltaDistY(ii,jj) = mean(deltaDistYTemp);
        
        %x change:
        deltaDataX(ii,jj) = deltaDataXTemp/numel(deltaDistXTemp);
        %y change:
        deltaDataY(ii,jj) = deltaDataYTemp/numel(deltaDistYTemp); 

        %x gradient:
        xGrad(ii,jj) = deltaDataXTemp/sum(deltaDistXTemp);
        %y gradient:
        yGrad(ii,jj) = deltaDataYTemp/sum(deltaDistYTemp);
        
        if isnan(xGrad(ii,jj))
            warning('gradientUtm:Nan',['gradient nan at indice: ' num2str(ii), ',' num2str(jj)])
        end
%         if isinf(xGrad(ii,jj))
%             keyboard
%             disp(['gradient nan at indice: ' num2str(ii), ',' num2str(jj)])
%         end
    end
end


if blPer == 1
    xGrad = xGrad./data;
    yGrad = yGrad./data;
    deltaDataX = deltaDataX./data;
    deltaDataY = deltaDataY./data;
end

%Set infinity values to nan:
xGrad(isinf(xGrad)) = nan;
yGrad(isinf(yGrad)) = nan;


%Assign output arguments
if nargout == 1 
    if blMag == 1
        varargout{1} = sqrt(xGrad.^2 + yGrad.^2);
    else
        error('gradientUtm:wrongOutZ', ['The number of output arguments '...
            'should be one because the z gradient option was selected.']);
    end
elseif nargout == 2 && blMag == 0
    varargout{1} = xGrad;
    varargout{2} = yGrad;
elseif nargout == 3
    varargout{1} = xGrad;
    varargout{2} = yGrad;
    if blMag == 1
        varargout{3} = sqrt(xGrad.^2 + yGrad.^2);
    else
        varargout{3} = sqrt(deltaDataX.^2 + deltaDataY.^2);
    end
elseif nargout == 4 && blMag == 0
    varargout{1} = xGrad;
    varargout{2} = yGrad;
    varargout{3} = deltaDataX;
    varargout{4} = deltaDataY;
elseif nargout == 5
    varargout{1} = xGrad;
    varargout{2} = yGrad;
    varargout{3} = deltaDataX;
    varargout{4} = deltaDataY;
    if blMag == 1
        varargout{5} = sqrt(xGrad.^2 + yGrad.^2);
    else
        varargout{5} = sqrt(deltaDataX.^2 + deltaDataY.^2);
    end
elseif nargout == 6 && blMag == 0  
    varargout{1} = xGrad;
    varargout{2} = yGrad;
    varargout{3} = deltaDataX;
    varargout{4} = deltaDataY;
    varargout{5} = deltaDistX;
    varargout{6} = deltaDistY;
elseif nargout == 7  
    varargout{1} = xGrad;
    varargout{2} = yGrad;
    varargout{3} = deltaDataX;
    varargout{4} = deltaDataY;
    varargout{5} = deltaDistX;
    varargout{6} = deltaDistY;
    if blMag == 1
        varargout{7} = sqrt(xGrad.^2 + yGrad.^2);
    else
        varargout{7} = sqrt(deltaDataX.^2 + deltaDataY.^2);
    end
else
    error('gradientUtm:wrongNOut', ['The number of output arguments '...
        'should be 1, 2, 4, or 6.']);
end
%For testing:
%[testX, testY] = gradient(data);
%testRatio = (testX/sum2d(testX))./(xGrad/sum2d(xGrad));