% STAGE 4b (simulation): OQPSK link with Doppler S-curve, timing offset, clock
% mismatch, optional tone interferer -> full OQPSK receiver. Validates the chain
% BEFORE it is pointed at a real Meteor-M recording.
% beta = 0.6 is the roll-off usually quoted for Meteor LRPT -- verify for your satellite.
clear; clc; close all; rng(4);
cfg = struct('sps',8,'beta',0.6,'span',10,'Rs',72e3, ...
             'fmax',3.2e3,'tauSym',0.3,'ppm',20,'supp','none');
%% 1) Single run: Doppler tracking + constellation
out = run_link_oqpsk(12, 1e5, cfg);
fprintf('Coarse CFO %.0f Hz | BER @12 dB: %.2e (%d bits)\n', out.f0, out.ber, out.nbits);
M = min(numel(out.fest), numel(out.ftrue)); t = (0:M-1)/cfg.Rs;
figure; plot(t, out.ftrue(1:M)/1e3, 'c', 'LineWidth', 2); hold on
plot(t, out.fest(1:M)/1e3, 'r--', 'LineWidth', 1.2); grid on
xlabel('Time (s)'); ylabel('kHz'); legend('True Doppler','Estimated (coarse + loop)'); title('OQPSK: Doppler tracking');
figure; r = max(1,numel(out.y)-1e4):numel(out.y);
plot(real(out.y(r)), imag(out.y(r)), '.'); axis square; grid on
title('OQPSK constellation after joint timing/carrier loop (I and Q sampled T/2 apart)');

%% 2) BER vs Eb/N0 (OQPSK has the same theoretical BER as QPSK)
EbN0 = 2:2:10; ber = zeros(size(EbN0)); nb = ber;
for k = 1:numel(EbN0), o = run_link_oqpsk(EbN0(k), 5e4, cfg); ber(k) = o.ber; nb(k) = o.nbits; end
figure; semilogy(EbN0, max(ber, 1./nb), 'o-', EbN0, 0.5*erfc(sqrt(10.^(EbN0/10))), 'k--'); grid on
xlabel('E_b/N_0 (dB)'); ylabel('BER'); legend('OQPSK full chain (Doppler+timing+clock)','Theory'); ylim([1e-5 1]);
title('OQPSK BER (floor = 1 error resolution)');

%% 3) With tone interferer (SIR = -10 dB) and the Stage 3 IIR notch
cfg.intf = struct('f',25e3,'sir_dB',-10);
for m = {'none','iir'}
    cfg.supp = m{1}; o = run_link_oqpsk(10, 5e4, cfg);
    fprintf('Tone interferer, suppression = %-4s : BER = %.2e\n', m{1}, o.ber);
end
