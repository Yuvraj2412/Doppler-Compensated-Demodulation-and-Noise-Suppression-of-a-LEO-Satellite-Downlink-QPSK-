% STAGE 4d-2: does the demodulated stream really contain Meteor LRPT frames?
% Takes ~1 s of symbols from the recording, runs the Stage 4c QPSK receiver, then
% Viterbi-decodes (r=1/2, K=7) under all 8 symbol mappings and searches for the
% 32-bit sync marker 0x1ACFFC1D repeating every 8192 bits. Needs meteor_bb.mat.
clear; clc; close all;
%% ---- CONFIG -------------------------------------------------------------
file     = 'C:\Users\tomar\Downloads\DSP_Project\LRPT_16-29-02_137900kHz.wav';
startSec = 30;          % where to take the segment (pick a strong-signal part)
beta     = 0.6;  loopBW = [0.003 0.005];  Nf = 2^16;
nSym     = 65536;       % symbols to decode (8 frames of 8192)
%% -------------------------------------------------------------------------
S = load('meteor_bb.mat','f_off','Rs','sps'); f_off = S.f_off; Rs = S.Rs; sps = S.sps;
info = audioinfo(file); Fs_in = info.SampleRate; Fs_out = Rs*sps;
durSec = (nSym + 8000)/Rs + 0.1;
x  = load_iq(file, 'wav', Fs_in, startSec, durSec);  x = x - mean(x);
xb = to_baseband(x, Fs_in, f_off, Fs_out);  xb = xb / sqrt(mean(abs(xb).^2));
r  = real_rx_chain(xb, Rs, sps, beta, 'qpsk', loopBW, Nf);
fprintf('Receiver: EVM = %.1f %%, %d symbols recovered\n', r.evm, numel(r.y));
y = r.y(4001:end);  y = y(1:min(numel(y), nSym));

tic; R = lrpt_sync_check(y); fprintf('Viterbi search done in %.0f s\n\n', toc);
fprintf('variant (swap,rot) | re-encode symbol errors | best sync match (of 32) | hits>=30 | repeating @8192\n');
for v = 1:numel(R)
    fprintf('      (%d,%d)        |        %5.1f %%          |          %2d             |   %3d    |   %3d\n', ...
            R(v).sw, R(v).rot, 100*R(v).err, R(v).best, R(v).hits, R(v).periodic);
end
[~,vb] = min([R.err]);
fprintf('\nBest variant: swap=%d, rot=%d -> re-encode symbol error %.2f %%\n', R(vb).sw, R(vb).rot, 100*R(vb).err);
if R(vb).periodic >= 2
    fprintf('SYNC MARKER FOUND repeating every 8192 bits (%d repeats): the stream is valid LRPT-style CCSDS framing.\n', R(vb).periodic);
else
    fprintf('No repeating sync marker found. See notes: segment too noisy, slips, or different code/frame convention.\n');
end
figure; subplot(2,1,1); bar(100*[R.err]); grid on; ylabel('%'); title('Re-encode symbol errors per variant (~12 % = wrong mapping, ~0 = right one)');
set(gca,'XTickLabel', arrayfun(@(q) sprintf('s%d r%d', q.sw, q.rot), R, 'UniformOutput', false));
subplot(2,1,2); plot(R(vb).M); grid on; xlabel('Decoded bit offset'); ylabel('Sync matches (of 32)');
title('Sync-marker correlation, best variant (peaks should repeat every 8192 bits)');
