function varargout = gradient_map(varargin)
%Exact same as built-in gradient function, but assumes that North is top of
%array and should be positive

if nargout == 2 
    if nargin == 1
        [tempX, tempY] = gradient(varargin{1});
    elseif nargin == 3
        [tempX, tempY] = gradient(varargin{1}, varargin{2}, varargin{3});
    else
        error('gradientMap:unknownInput','This function is only designed to work with 1 or 3 input arguments.');
    end
    
    varargout = cell(2,1);
    varargout{1} = tempX;
    varargout{2} = tempY;
    if nargin == 1 %Only necessary to flip when coordinates are not provided.
        varargout{2} = -varargout{2};
    end
else
    error('gradientMap:unknownOutput','This function is only designed to work with two output arguments.');
end