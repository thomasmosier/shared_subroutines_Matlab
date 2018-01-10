function y = cell_mode(x)

y = unique(x);
n = zeros(length(y), 1);
for iy = 1:length(y)
  n(iy) = length(find(strcmp(y{iy}, x)));
end
[~, itemp] = max(n);
y = y(itemp);