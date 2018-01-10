function q = alignment_factor(grid1x, grid1y, grid2x, grid2y)

[th1, ~] = cart2pol(grid1x, grid1y);
[th2, ~] = cart2pol(grid2x, grid2y);
%Calculate cosine of relative alignment:
q = cos(th1-th2);