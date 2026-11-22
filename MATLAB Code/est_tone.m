function [fi, snr_dB] = est_tone(x, Fs)
% Strongest narrowband component: frequency (Hz) and detection strength (dB above
% the local spectral level). Pure noise/signal gives ~10-14 dB; a real tone gives far more.
N = min(numel(x), 2^17); NF = 2^18;
w = 0.5 - 0.5*cos(2*pi*(0:N-1).'/N);
P = abs(fft(x(1:N).*w, NF)).^2 + eps;  Lg = log(P);
[~,k] = max(P);
a = Lg(mod(k-2,NF)+1); b = Lg(k); c = Lg(mod(k,NF)+1);
d = 0.5*(a - c)/(a - 2*b + c);
fi = (k-1+d)/NF*Fs;  if fi >= Fs/2, fi = fi - Fs; end
idx = mod([k-300:k-20, k+20:k+300] - 1, NF) + 1;
snr_dB = 10*log10(P(k)/mean(P(idx)));
end
