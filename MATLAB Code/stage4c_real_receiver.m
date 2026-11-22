% STAGE 4c: run the receiver on the baseband file made by Stage 4a (meteor_bb.mat).
% No transmitted bits exist, so quality is judged by EVM and the constellation.
% It tries QPSK and OQPSK, and two RRC roll-offs, and keeps the lowest-EVM result.
clear; clc; close all;
load('meteor_bb.mat');                          % xb, Fs_out, Rs, sps
useSec  = 4;                                     % seconds of signal to analyse
startSec = 0;                                    % where in meteor_bb to start
loopBW  = [0.003 0.005];                         % [timing carrier]; widen if Doppler loses lock
Nf      = 2^16;                                  % coarse-CFO FFT length (samples)
modes = {'qpsk','oqpsk'};  betas = [0.35 0.6];

a = round(startSec*Fs_out)+1; b = min(numel(xb), a + round(useSec*Fs_out) - 1);
x = xb(a:b);  fprintf('Analysing %.2f s (%d samples) at Rs = %.0f sym/s\n', numel(x)/Fs_out, numel(x), Rs);

best = []; k = 0; R = {};
for i = 1:numel(modes)
    for j = 1:numel(betas)
        k = k + 1; R{k} = real_rx_chain(x, Rs, sps, betas(j), modes{i}, loopBW, Nf);
        fprintf('%-5s beta=%.2f : EVM = %5.1f %%   CFO = %7.1f Hz   4th-power line = %4.1f dB\n', ...
                modes{i}, betas(j), R{k}.evm, R{k}.f0, R{k}.line_dB);
        if isempty(best) || R{k}.evm < best.evm, best = R{k}; end
    end
end
fprintf('\nBest: %s, beta = %.2f, EVM = %.1f %%\n', upper(best.mode), best.beta, best.evm);
fprintf('(For reference: clean 4-cluster QPSK at ~10 dB Es/N0 is roughly 25-30 %% EVM; random noise ~ 70+ %%.)\n');

figure; t = (0:numel(best.fest)-1)/Rs + startSec;
plot(t, best.fest/1e3); grid on; xlabel('Time (s)'); ylabel('Frequency offset (kHz)');
title(sprintf('Estimated carrier offset (%s)', upper(best.mode)));
figure; r = max(1, numel(best.y)-2e4):numel(best.y);
plot(real(best.y(r)), imag(best.y(r)), '.'); axis square; grid on; xlim([-2 2]); ylim([-2 2]);
title(sprintf('Recovered constellation: %s, EVM %.1f%%', upper(best.mode), best.evm));
figure; bar(cellfun(@(q) q.evm, R)); grid on; ylabel('EVM (%)');
set(gca,'XTickLabel', arrayfun(@(q) sprintf('%s %.2f', q{1}.mode, q{1}.beta), R, 'UniformOutput', false));
title('Receiver configuration comparison (lower is better)');
