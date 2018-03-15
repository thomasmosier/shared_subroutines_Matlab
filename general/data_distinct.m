function [blDist, pVal] = data_distinct(data1, data2, nBt, sig)

%Using bootstrapping loop to test significance of surrogates
matKStat = nan([nBt, 1], 'single');

%NOTE: Andy mentioned alternative formulation using the Lilliefors test formulation (w/o FFT scrambling?).
%Loop over bootstrapping iterations
for mm = 1 : nBt
   %(1) Take the FFT of the data, and randomize the phase of each component: 
   temp1F = fft(data1); 
   temp1FPr = temp1F.*exp(2*pi*1i*rand(size(temp1F)));

   temp2F = fft(data2); 
   temp2FPr = temp2F.*exp(2*pi*1i*rand(size(temp2F)));

   %(2) Transform output of #1 back to the time domain. This gives a surrogate time series:
   temp1FPr0 = ifft(temp1FPr,'symmetric');
   temp1Pr = (temp1FPr0-nanmean(temp1FPr0))/nanstd(temp1FPr0);

   temp2FPr0 = ifft(temp2FPr,'symmetric');
   temp2Pr = (temp2FPr0-nanmean(temp2FPr0))/nanstd(temp2FPr0);

   %(3) Determine significance using the surrogate:
   %Returns 1 if the distributions are distinct, 0 otherwise
   [~,~, matKStat(mm)] = kstest2(temp1Pr, temp2Pr);
end

%(5) Now, you have a distribution of the test statistic that takes into account autocorrelation. The adjusted critical value is just:
kstatCrit = prctile(matKStat, round(100*(1-sig)));
% If you now perform the test on the original time series...
[~,~, kstatY] = kstest2(data1, data2);
% ... reject the null if the result is greater than or equal to this critical value:
blDist = kstatY >= kstatCrit;

% P-value estimate (this gets more accurate as nboot increases):
[F,x] = ecdf(matKStat);
[~,idx] = min(abs(kstatY - x));
pVal = 1-F(idx);