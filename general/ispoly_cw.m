function tf = ispoly_cw(xCl,yCl)

if ~isequal(class(xCl), class(yCl))
    error('ispolyCw:classMismatch',['The x input is class ' class(xCl) ' and the y input is class ' class(yCl) '. These must be the same.']);
end

%Check if polygons need to be split
if isnumeric(xCl) && any(isnan(xCl))
    [xCl, yCl] = poly_split(xCl, yCl);
end

if isnumeric(xCl)
    xCl = {xCl};
    yCl = {yCl};
    
end

tf = nan(numel(xCl),1);

for ii = 1 : numel(xCl(:))
    if any(isnan(xCl{ii})) || any(isnan(yCl{ii}))
        warning('ispolyCw:cellNan', ['Polygon ' num2str(ii) ' contains nan values. This is not allowed.']);
        continue
    end
    
    if numel(xCl{ii}) <= 2
        tf(ii) = 1;
    end
    
    
    %Calculate signed area (https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order)
    %A = 1/2 * (x1*y2 - x2*y1 + x2*y3 - x3*y2 + ... + xn*y1 - x1*yn)
    a = 0;
    for jj = 1 : numel(xCl{ii}) - 1
        a = a + xCl{ii}(jj)*yCl{ii}(jj+1) - xCl{ii}(jj+1)*yCl{ii}(jj);
    end
    clear jj
    a = a + xCl{ii}(end)*yCl{ii}(1) - xCl{ii}(1)*yCl{ii}(end);
    
    if a > 0
        tf(ii) = 0;
    elseif a < 0
        tf(ii) = 1;
    end
    
        
%     %Start by finding lowest y with smallest x:
%     [~, indLl] = min(yCl{ii});
%     indLl = find(yCl{ii} == yCl{ii}(indLl));
%     if numel(indLl) > 1
%         [~, indXTemp] = min(xCl{ii}(indLl));
%         indXTemp = find(xCl{ii} == xCl{ii}(indXTemp));
%         if numel(indXTemp) > 1 %This is indeterminant (for this simple version of function)
%             continue
%         end
%         indLl = indLl(indXTemp);
%     end
end