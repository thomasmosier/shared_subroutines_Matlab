function [f1, P1, f2, P2] = fft_amp(x, dT, varargin)

nrm = 0;
pwrSpec = 0;
for ii = 1 : numel(varargin(:))
    if strcmpi(varargin{ii}, 'norm') || strcmpi(varargin{ii}, 'normalize')
       nrm = 1; 
    elseif strcmpi(varargin{ii}, 'power') 
        pwrSpec = 1;
    end
end

Fs = 1/dT; %Sampling frequency

L = numel(x);
f1 = (0:(L/2))'*(Fs/L);
f2 = (0: L-1 )'*(Fs/L);  


pwr = fft(x);

if pwrSpec == 1
    pwr = pwr.^2;
end

P2 = abs(pwr/L);
P1 = P2(int16(1:L/2+1));
P1(2:end-1) = 2*P1(2:end-1);

if nrm == 1
    P1 = P1/max(P1);
    P2 = P2/max(P2);
%    pwr = fftshift(pwr/numel(pwr)); 
end