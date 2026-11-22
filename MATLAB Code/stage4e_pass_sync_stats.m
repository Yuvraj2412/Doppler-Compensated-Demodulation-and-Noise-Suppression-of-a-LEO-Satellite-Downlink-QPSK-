% STAGE 4e: frame-sync statistics over the whole recording.
% For each chunk: QPSK receiver -> Viterbi under all 8 mappings -> sync-marker search.
% Reports, per chunk, the raw symbol error rate (re-encode check) and the number of
% 1024-byte frames whose 0x1ACFFC1D marker was found. Needs meteor_bb.mat (Stage 4a).
clear; clc; close all;
%% ---- CONFIG -------------------------------------------------------------
file     = 'C:\Users\tomar\Downloads\DSP_Project\LRPT_16-29-02_137900kHz.wav';
chunkSec = 1.4;          % signal per chunk (gives >65536 symbols after loop settling)
hopSec   = 3;            % spacing between chunk starts (use 1.4 for full coverage)
beta     = 0.6;  loopBW = [0.003 0.005];  Nf = 2^16;
nDec     = 65536;        % symbols decoded per chunk (8 frames)
syncMin  = 6;            % chunk counts as frame-synchronised if >= this many markers found
%% -------------------------------------------------------------------------
S = load('meteor_bb.mat','f_off','Rs','sps'); f_off = S.f_off; Rs = S.Rs; sps = S.sps;
info = audioinfo(file); Fs_in = info.SampleRate; Ttot = info.Duration; Fs_out = Rs*sps;
starts = 0:hopSec:(Ttot - chunkSec);  nC = numel(starts);  tc = starts + chunkSec/2;
evm = nan(1,nC); err = nan(1,nC); hits = zeros(1,nC); sw = nan(1,nC); rot = nan(1,nC); inv = nan(1,nC);
fprintf('%d chunks of %.1f s from a %.1f s recording\n', nC, chunkSec, Ttot); tic
for i = 1:nC
    x  = load_iq(file, 'wav', Fs_in, starts(i), chunkSec);  x = x - mean(x);
    xb = to_baseband(x, Fs_in, f_off, Fs_out);  xb = xb / sqrt(mean(abs(xb).^2));
    try
        r = real_rx_chain(xb, Rs, sps, beta, 'qpsk', loopBW, Nf);  evm(i) = r.evm;
        y = r.y(4001:end);
        if numel(y) < 20000, error('too few symbols'); end
        y = y(1:min(numel(y), nDec));
        R = lrpt_sync_check(y);  e = [R.err];
        good = find(e <= min(e) + 0.002);                       % ties = 180-degree pair
        pick = good(find(~[R(good).inv], 1)); if isempty(pick), pick = good(1); end
        err(i) = R(pick).err; hits(i) = R(pick).hits; sw(i) = R(pick).sw; rot(i) = R(pick).rot; inv(i) = R(pick).inv;
    catch ex
        fprintf('  chunk %d failed: %s\n', i, ex.message);
    end
    fprintf('chunk %2d/%d t=%6.1f s  EVM %5.1f %%  mapping (%d,%d)  raw sym err %5.2f %%  markers found: %d   (%.0f s)\n', ...
            i, nC, tc(i), evm(i), sw(i), rot(i), 100*err(i), hits(i), toc);
end
sync = hits >= syncMin;
fprintf('\nFrame-synchronised chunks: %d of %d\n', sum(sync), nC);
fprintf('Markers found: %d total, i.e. %.1f per chunk (8 expected for %d decoded symbols)\n', sum(hits), mean(hits), nDec);
fprintf('Median raw symbol error (synchronised chunks): %.3f %%\n', 100*median(err(sync)));
fprintf('Mappings seen (swap,rot): '); u = unique([sw(sync); rot(sync)]','rows'); disp(u');

figure;
subplot(3,1,1); plot(tc, hits, 'o-'); grid on; ylabel('Markers found'); title('Sync markers per chunk (8 = every frame found)');
subplot(3,1,2); semilogy(tc, max(100*err, 1e-3), 'o-'); grid on; ylabel('Raw symbol error (%)'); title('Re-encode symbol error per chunk (floor 0.001 %)');
subplot(3,1,3); plot(tc, evm, 'o-'); grid on; ylabel('EVM (%)'); xlabel('Time (s)'); title('EVM per chunk');
save('pass_sync_stats.mat','tc','evm','err','hits','sw','rot','inv');
