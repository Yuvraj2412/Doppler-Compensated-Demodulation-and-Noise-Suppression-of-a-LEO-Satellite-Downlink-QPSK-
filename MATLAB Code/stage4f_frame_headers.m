% STAGE 4f: de-scramble the decoded frames and read the VCDU headers.
% If everything upstream is right, the 24-bit frame counter must increase by exactly 1
% from frame to frame, and by (frame rate x elapsed time) across the whole pass
% (frame rate = 72000 sym/s / 8192 = 8.789 frames/s). Needs meteor_bb.mat (Stage 4a).
clear; clc; close all;
%% ---- CONFIG -------------------------------------------------------------
file     = 'C:\Users\tomar\Downloads\DSP_Project\LRPT_16-29-02_137900kHz.wav';
chunkSec = 1.4;  hopSec = 3;
beta     = 0.6;  loopBW = [0.003 0.005];  Nf = 2^16;  nDec = 65536;
%% -------------------------------------------------------------------------
pn = ccsds_pn_bytes(1020);
assert(isequal(pn(1:4), [255 72 14 192]), 'PN generator check failed');
asm = double(dec2bin(hex2dec('1ACFFC1D'), 32).' == '1');
S = load('meteor_bb.mat','f_off','Rs','sps'); f_off = S.f_off; Rs = S.Rs; sps = S.sps;
info = audioinfo(file); Fs_in = info.SampleRate; Ttot = info.Duration; Fs_out = Rs*sps;
starts = 0:hopSec:(Ttot - chunkSec);  nC = numel(starts);
T = []; CNT = []; SC = []; VC = []; SG = []; CH = [];         % one row per frame
fprintf('%d chunks\n', nC); tic
for i = 1:nC
    x  = load_iq(file, 'wav', Fs_in, starts(i), chunkSec);  x = x - mean(x);
    xb = to_baseband(x, Fs_in, f_off, Fs_out);  xb = xb / sqrt(mean(abs(xb).^2));
    try
        r = real_rx_chain(xb, Rs, sps, beta, 'qpsk', loopBW, Nf);
        y = r.y(4001:end);  y = y(1:min(numel(y), nDec));
        R = lrpt_sync_check(y);  e = [R.err];  good = find(e <= min(e) + 0.002);
        pick = good(find(~[R(good).inv], 1)); if isempty(pick), pick = good(1); end
        u = double(R(pick).u);  if R(pick).inv, u = 1 - u; end
        Md = zeros(numel(u)-31, 1);  for k = 1:32, Md = Md + (u(k:end-32+k) == asm(k)); end
        pos = find(Md >= 30).';
        for p = pos
            if p + 8191 > numel(u), continue; end
            by = reshape(u(p:p+8191), 8, []).' * 2.^(7:-1:0).';  by = by(:).';
            d  = bitxor(by(5:1024), pn);
            T(end+1) = starts(i) + (4000 + p)/Rs;                    % approx. time of frame start
            SC(end+1) = bitor(bitshift(bitand(d(1),63),2), bitshift(d(2),-6));
            VC(end+1) = bitand(d(2),63);  CNT(end+1) = d(3)*65536 + d(4)*256 + d(5);
            SG(end+1) = d(6);  CH(end+1) = i;
        end
    catch ex
        fprintf('  chunk %d failed: %s\n', i, ex.message);
    end
    fprintf('chunk %2d/%d  frames so far: %d   (%.0f s)\n', i, nC, numel(CNT), toc);
end
fprintf('\nFrames read: %d\n', numel(CNT));
[u1,~,j1] = unique([SC; VC]', 'rows'); n1 = accumarray(j1, 1);
for k = 1:size(u1,1), fprintf('  SCID %3d, VCID %2d : %d frames\n', u1(k,1), u1(k,2), n1(k)); end
same = CH(2:end) == CH(1:end-1);  dc = diff(CNT);
fprintf('Consecutive frames inside a chunk: %d pairs, counter step = +1 in %.1f %%\n', sum(same), 100*mean(dc(same) == 1));
if numel(CNT) > 10
    ok = (SC == mode(SC)) & (VC == mode(VC));
    p = polyfit(T(ok), CNT(ok), 1);
    fprintf('Counter vs time: %.4f frames/s (nominal 72000/8192 = %.4f), residual rms = %.2f frames\n', p(1), 72000/8192, std(CNT(ok) - polyval(p, T(ok))));
    figure; subplot(2,1,1); plot(T(ok), CNT(ok), '.'); grid on; ylabel('Frame counter'); title('VCDU frame counter vs time (straight line = no lost frames)');
    subplot(2,1,2); plot(T(ok), CNT(ok) - polyval(p, T(ok)), '.'); grid on; xlabel('Time (s)'); ylabel('Residual (frames)'); title('Residual after linear fit');
end
save('pass_headers.mat','T','CNT','SC','VC','SG','CH');
