% STAGE 1: clean QPSK link + AWGN, RRC matched filter, BER vs theory.
clear; clc; close all; rng(1);
sps = 8; beta = 0.35; span = 10; Rs = 72e3; Fs = Rs*sps; L = span*sps+1;
h = rrc_filter(beta, span, sps);
EbN0 = 0:1:10; nSym = 2e5; ber = zeros(size(EbN0));
for k = 1:numel(EbN0)
    bits = randi([0 1], 2*nSym, 1);  s = qpsk_mod(bits);
    up = zeros(nSym*sps,1); up(1:sps:end) = s;
    tx = filter(h,1,[up; zeros(L-1,1)]);
    sigma2 = 0.5/10^(EbN0(k)/10);                    % Es = 1
    rx = tx + sqrt(sigma2/2)*(randn(size(tx)) + 1j*randn(size(tx)));
    y  = filter(h,1,[rx; zeros(L-1,1)]);
    y  = y(L : sps : L+(nSym-1)*sps);                % optimum sampling instants
    ber(k) = mean(qpsk_demod(y) ~= bits);
end
theory = 0.5*erfc(sqrt(10.^(EbN0/10)));
figure; semilogy(EbN0, ber, 'o-', EbN0, theory, 'k--', 'LineWidth', 1.2); grid on
xlabel('E_b/N_0 (dB)'); ylabel('BER'); legend('Simulated','Theory Q(\surd(2E_b/N_0))');
title('Stage 1: QPSK over AWGN with RRC pulse shaping'); ylim([1e-6 1]);

% Spectrum + eye diagram sanity checks
[P,f] = pwelch_simple(tx, 1024, Fs);
figure; subplot(1,2,1); plot(f/1e3, 10*log10(P)); grid on; xlabel('kHz'); ylabel('dB'); title('TX PSD');
subplot(1,2,2); e = real(tx(L:L+2000*sps)); e = reshape(e(1:floor(numel(e)/(2*sps))*2*sps), 2*sps, []);
plot(0:2*sps-1, e(:,1:300), 'b'); grid on; title('Eye diagram (I)');

function [P,f] = pwelch_simple(x, N, Fs)   % local fn: needs MATLAB R2016b+
w = 0.5 - 0.5*cos(2*pi*(0:N-1).'/N);  K = floor(numel(x)/N); P = zeros(N,1);
for i = 1:K, P = P + abs(fftshift(fft(x((i-1)*N+1:i*N).*w))).^2; end
P = P/(K*sum(w.^2)*Fs); f = (-N/2:N/2-1).'*Fs/N;
end
