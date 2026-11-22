% STAGE 4d-1: whole-recording analysis. Processes the WAV in short chunks and
% tracks, over the full pass: carrier frequency (Doppler), EVM and the EVM-implied
% SNR. Needs meteor_bb.mat from Stage 4a (for f_off, Rs, sps) in the current folder.
clear; clc; close all;
%% ---- CONFIG -------------------------------------------------------------
file     = 'C:\Users\tomar\Downloads\DSP_Project\LRPT_16-29-02_137900kHz.wav';
chunkSec = 2;            % length of each analysed chunk (s)
hopSec   = 4;            % spacing between chunk starts (s); = chunkSec for 100% coverage
beta     = 0.6;          % RRC roll-off (0.35 and 0.6 were indistinguishable in 4c)
loopBW   = [0.003 0.005];% [timing carrier] loop BnT
Nf       = 2^16;         % coarse-CFO FFT length
evmLock  = 40;           % chunks with EVM below this (%) count as "locked"
%% -------------------------------------------------------------------------
if exist('meteor_bb.mat','file'), S = load('meteor_bb.mat','f_off','Rs','sps');
    f_off = S.f_off; Rs = S.Rs; sps = S.sps;
else, f_off = 0; Rs = 72e3; sps = 8; warning('meteor_bb.mat not found: using f_off = 0'); end
info = audioinfo(file); Fs_in = info.SampleRate; Ttot = info.Duration; Fs_out = Rs*sps;
fprintf('File: %.1f s at %.0f Sa/s, %d channels. f_off = %.1f Hz\n', Ttot, Fs_in, info.NumChannels, f_off);
starts = 0:hopSec:(Ttot - chunkSec);  nC = numel(starts);
tc = starts + chunkSec/2;  evm = nan(1,nC); fabs = nan(1,nC); line_dB = nan(1,nC);
tic
for i = 1:nC
    x  = load_iq(file, 'wav', Fs_in, starts(i), chunkSec);  x = x - mean(x);
    xb = to_baseband(x, Fs_in, f_off, Fs_out);  xb = xb / sqrt(mean(abs(xb).^2));
    try
        r = real_rx_chain(xb, Rs, sps, beta, 'qpsk', loopBW, Nf);
        evm(i) = r.evm;  line_dB(i) = r.line_dB;
        fabs(i) = f_off + mean(r.fest(round(end/5):end));      % carrier in file frame (Hz)
    catch err
        fprintf('  chunk %d failed: %s\n', i, err.message);
    end
    fprintf('chunk %2d/%d  t=%6.1f s  EVM = %5.1f %%  carrier = %8.1f Hz   (%.0f s elapsed)\n', ...
            i, nC, tc(i), evm(i), fabs(i), toc);
end
snr = -20*log10(evm/100);                      % EVM-implied SNR (rough, biased at low SNR)
ok = evm < evmLock;
fprintf('\nLocked chunks: %d of %d (EVM < %d %%). Median EVM (locked) = %.1f %%\n', sum(ok), nC, evmLock, median(evm(ok)));

figure;
subplot(3,1,1); plot(tc(ok), fabs(ok)/1e3, 'o-'); hold on; grid on
if any(~ok), plot(tc(~ok), fabs(~ok)/1e3, 'rx'); end
ylabel('Carrier (kHz)'); title('Carrier frequency in file frame (Doppler + fixed offsets)');
if sum(ok) >= 4
    p = polyfit(tc(ok), fabs(ok), 3); fit = polyval(p, tc);
    plot(tc, fit/1e3, 'k--'); legend('Locked chunks','Cubic fit','Location','northeast'); slope = polyval(polyder(p), tc);
    fprintf('Cubic fit: change %.0f Hz over %.0f s, max |slope| = %.1f Hz/s, residual rms = %.1f Hz\n', ...
            fit(end)-fit(1), tc(end)-tc(1), max(abs(slope)), std(fabs(ok) - polyval(p, tc(ok))));
end
subplot(3,1,2); plot(tc, evm, 'o-'); hold on; plot([tc(1) tc(end)], [evmLock evmLock], 'r--'); grid on; ylabel('EVM (%)');
title('EVM per chunk (lower is better)');
subplot(3,1,3); plot(tc, snr, 'o-'); grid on; ylabel('SNR (dB)'); xlabel('Time (s)'); title('EVM-implied SNR');
save('pass_results.mat','tc','evm','fabs','snr','line_dB','Rs','f_off');
