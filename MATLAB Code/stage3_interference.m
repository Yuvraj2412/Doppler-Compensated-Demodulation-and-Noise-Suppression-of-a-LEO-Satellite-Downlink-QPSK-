% STAGE 3: narrowband interferer + suppression (IIR notch, FIR band-stop,
% NLMS tone canceller, NLMS adaptive line enhancer) on top of Stage 2 chain.
% The interferer is a fixed-frequency tone at +25 kHz (inside the signal band),
% SIR = -10 dB. It is NOT Doppler shifted (ground-based source), while the
% satellite signal is.
clear; clc; close all; rng(3);
cfg = struct('sps',8,'beta',0.35,'span',10,'Rs',72e3, ...
             'fmax',3.2e3,'tauSym',0.3,'ppm',20);
Fs = cfg.Rs*cfg.sps;  nSym = 5e4;
cfg.intf = struct('f',25e3,'sir_dB',-10);
modes  = {'none','iir','fir','nlms','ale'};
labels = {'No suppression','IIR notch','FIR band-stop (511 tap)','NLMS tone canceller','NLMS line enhancer'};

%% A) Spectra before / after suppression (Eb/N0 = 10 dB)
res = cell(1,numel(modes));
for i = 1:numel(modes)
    cfg.supp = modes{i};  res{i} = run_link(10, nSym, cfg);
    fprintf('%-24s BER = %.2e   (est. tone %.1f Hz)\n', labels{i}, res{i}.ber, res{i}.fi);
end
figure; hold on; col = lines(5);
[S,f] = simple_spectrogram(res{1}.rx_pre, Fs, 2048, 1024);
plot(f/1e3, 10*log10(mean(10.^(S/10),2)), 'Color',[.6 .6 .6], 'LineWidth',1.5);
leg = {'Before (with interferer)'};
for i = 2:numel(modes)
    [S,f] = simple_spectrogram(res{i}.rx_post, Fs, 2048, 1024);
    plot(f/1e3, 10*log10(mean(10.^(S/10),2)), 'Color',col(i,:)); leg{end+1} = labels{i};
end
grid on; xlim([-100 100]); xlabel('Frequency (kHz)'); ylabel('PSD (dB/Hz)');
title('Received spectrum before and after interference suppression'); legend(leg,'Location','northeast');

%% B) Spectrograms: before vs after NLMS tone canceller
figure;
[S,f,t] = simple_spectrogram(res{1}.rx_pre,  Fs, 1024, 512);
subplot(1,2,1); imagesc(t, f/1e3, S); axis xy; ylim([-100 100]); colorbar
xlabel('Time (s)'); ylabel('kHz'); title('Spectrogram: before'); cl = caxis;
[S,f,t] = simple_spectrogram(res{4}.rx_post, Fs, 1024, 512);
subplot(1,2,2); imagesc(t, f/1e3, S); axis xy; ylim([-100 100]); caxis(cl); colorbar
xlabel('Time (s)'); ylabel('kHz'); title('Spectrogram: after NLMS canceller');

%% C) BER vs Eb/N0 at SIR = -10 dB
EbN0 = 2:2:10;  ber = zeros(numel(modes), numel(EbN0));  nb = ber;
for i = 1:numel(modes)
    cfg.supp = modes{i};
    for k = 1:numel(EbN0)
        o = run_link(EbN0(k), nSym, cfg); ber(i,k) = o.ber; nb(i,k) = o.nbits;
    end
end
figure; semilogy(EbN0, 0.5*erfc(sqrt(10.^(EbN0/10))), 'k--', 'LineWidth',1.5); hold on
for i = 1:numel(modes), semilogy(EbN0, max(ber(i,:), 1./nb(i,:)), 'o-', 'Color', col(i,:)); end
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER'); ylim([1e-5 1]);
legend([{'Ideal AWGN (no interferer)'}, labels], 'Location','southwest');
title('BER with tone interferer, SIR = -10 dB (floor = 1 error resolution)');

%% D) BER vs SIR at Eb/N0 = 10 dB, plus cost of each filter with NO interferer
SIR = [-20 -10 0 10];  berS = zeros(numel(modes), numel(SIR));
for i = 1:numel(modes)
    cfg.supp = modes{i};
    for k = 1:numel(SIR)
        cfg.intf.sir_dB = SIR(k);  o = run_link(10, nSym, cfg); berS(i,k) = o.ber;
    end
end
figure; for i = 1:numel(modes), semilogy(SIR, max(berS(i,:),1e-5), 'o-', 'Color', col(i,:)); hold on; end
grid on; xlabel('SIR (dB)'); ylabel('BER'); legend(labels,'Location','southwest'); ylim([1e-5 1]);
title('BER vs interference level (E_b/N_0 = 10 dB)');

cfg.intf = [];                                     % no interferer present
fprintf('\nCost of each filter with NO interferer (Eb/N0 = 10 dB):\n');
for i = 1:numel(modes)
    cfg.supp = modes{i};  o = run_link(10, nSym, cfg);
    fprintf('%-24s BER = %.2e\n', labels{i}, o.ber);
end
