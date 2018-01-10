function vecOut = rand_weighted(vecNumbs, vecWgts, N)

sumWgts = cumsum(vecWgts);

%Wgts should sum to 1
if round2(sumWgts(end), 3) ~= 1
    error('randweighted:vecWgtSum',['The weights must sum to 1 (currently the sum is ' num2str(sumWgts(end)) ').'])
end

vecOut = nan(N, 1);

randSelect = rand(N,1);
for ii = 1 : N
    vecOut(ii) = vecNumbs(find(randSelect(ii) <= sumWgts, 1, 'first'));
end

if any(isnan(vecOut))
    error('randWeighted:nanOut', [num2str(sum(isnan(vecOut))) ' output indices are nan. 0 are expected to be nan.']);
end