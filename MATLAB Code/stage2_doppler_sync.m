% STAGE 2: AWGN + Doppler S-curve + timing/clock offset, with full sync chain.
% NOTE: a real pass is ~10 min; here the S-curve is compressed into ~3 s of
% signal (so ~200x faster Doppler rate than reality) to keep run time small.
clear; clc; close all; rng(2);
cfg = struct('sps',8,'beta',0.35,'span',10,'Rs',72e3, ...
             'fmax',3.2e3,'tauSym',0.3,'ppm',20);
out = run_link(12, 2e5, cfg);
fprintf('Coarse CFO estimate: %.1f Hz | BER @12 dB: %.2e\n', out.f0, out.ber);

figure; M = min(numel(out.fest), numel(out.ftrue)); t = (0:M-1)/cfg.Rs;
plot(t, out.ftrue(1:M)/1e3, 'c', 'LineWidth', 2); hold on
plot(t, out.fest(1:M)/1e3, 'r--', 'LineWidth', 1.2);xlabel('Time (s)'); ylabel('kHz'); legend('True Doppler','Estimated (coarse + Costas)');
title('Doppler tracking');

figure; r = max(1,numel(out.y)-1e4):numel(out.y);
subplot(1,2,1); plot(real(out.y_pre(r)), imag(out.y_pre(r)), '.'); axis square; grid on; title('After timing only');
subplot(1,2,2); plot(real(out.y(r)),     imag(out.y(r)),     '.'); axis square; grid on; title('After Costas');

EbN0 = 4:2:12; ber = zeros(size(EbN0));
for k = 1:numel(EbN0), o = run_link(EbN0(k), 1e5, cfg); ber(k) = o.ber; end
figure; semilogy(EbN0, max(ber,1e-6), 'o-', EbN0, 0.5*erfc(sqrt(10.^(EbN0/10))), 'k--'); grid on
xlabel('E_b/N_0 (dB)'); ylabel('BER'); legend('Full chain (Doppler+timing)','Ideal AWGN');
