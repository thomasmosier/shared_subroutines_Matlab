function rank = vec_rank(vec)

rank = nan(size(vec));

for ii = 1 : numel(vec)
    rank(ii) = sum(vec(:) <= vec(ii)) / numel(vec);
end