function [S, f, t] = simple_spectrogram(x, Fs, nfft, hop)
% STFT power (dB), Hann window, no toolbox needed.
x = x(:); w = 0.5 - 0.5*cos(2*pi*(0:nfft-1).'/nfft);
K = floor((numel(x)-nfft)/hop) + 1; S = zeros(nfft, K);
for i = 1:K
    seg = x((i-1)*hop+1 : (i-1)*hop+nfft) .* w;
    S(:,i) = fftshift(abs(fft(seg)).^2);
end
S = 10*log10(S/(sum(w.^2)*Fs) + eps);
f = (-nfft/2:nfft/2-1).'*Fs/nfft; t = ((0:K-1)*hop + nfft/2)/Fs;
end
