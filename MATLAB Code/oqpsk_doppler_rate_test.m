% Reproduces Table "OQPSK receiver: effect of Doppler rate and loop bandwidth" in the report.
% Compares a fast (+-3.2 kHz) and slow (+-300 Hz) Doppler swing with three loop bandwidths.
% Takes several minutes. Results use ~1e5 symbols per entry.
clear; clc;
base = struct('sps',8,'beta',0.6,'span',10,'Rs',72e3,'tauSym',0.3,'ppm',20);
bws = [0.005 0.01; 0.002 0.004; 0.001 0.002];
for fm = [3.2e3 300]
    for k = 1:size(bws,1)
        cfg = base; cfg.fmax = fm; cfg.loopBW = bws(k,:);
        rng(21); o6 = run_link_oqpsk(6, 1e5, cfg);
        rng(22); o2 = run_link_oqpsk(2, 1e5, cfg);
        fprintf('fmax=%4.0f Hz  loops=[%.3f %.3f]:  BER@6dB=%.2e (theory 2.4e-3)   BER@2dB=%.2e (theory 3.8e-2)\n', ...
                fm, bws(k,1), bws(k,2), o6.ber, o2.ber);
    end
end
