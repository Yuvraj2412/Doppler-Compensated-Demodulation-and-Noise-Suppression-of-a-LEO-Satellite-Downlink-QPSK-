function r = real_rx_chain(xb, Rs, sps, beta, mode, loopBW, Nf)
% Blind receiver for a baseband recording (no known bits).
% mode: 'qpsk' (Gardner + Costas) or 'oqpsk' (joint loop). loopBW = [timing carrier] BnT.
Fs = Rs*sps; L = 10*sps + 1; xb = xb(:);
h = rrc_filter(beta, 10, sps);
m = filter(h, 1, [xb; zeros(L-1,1)]);
% coarse CFO: 4th-power line (QPSK and OQPSK both show it)
NF = 2^nextpow2(8*Nf); a0 = 8193; seg = m(a0:a0+Nf-1);
X = abs(fft(seg.^4, NF)); [pk,k] = max(X);
f4 = (k-1)/NF*Fs; if f4 > Fs/2, f4 = f4 - Fs; end
f0 = f4/4;  r.line_dB = 20*log10(pk/median(X));
m = m .* exp(-1j*2*pi*f0*(0:numel(m)-1).'/Fs);
% two-pass amplitude normalisation (loops assume ~unit-energy symbols)
sc = 1/sqrt(mean(abs(m).^2)*sps);
for pass = 1:2
    if strcmp(mode,'qpsk')
        ys = gardner_sync(m*sc, sps, loopBW(1), 0.707);
        [yc, fHz] = costas_qpsk(ys, loopBW(2), 0.707, Rs);
    else
        [yc, fHz] = oqpsk_sync(m*sc, sps, loopBW(1), loopBW(2), 0.707, Rs);
    end
    if pass == 1, tail = yc(round(end/2):end); sc = sc/sqrt(mean(abs(tail).^2)); end
end
skip = min(3000, floor(numel(yc)/4)); z = yc(skip+1:end);
ideal = (sign(real(z)) + 1j*sign(imag(z))) * mean(abs(real(z)) + abs(imag(z)))/2;
r.evm = 100*sqrt(mean(abs(z-ideal).^2)/mean(abs(ideal).^2));
r.y = yc; r.fest = f0 + fHz; r.f0 = f0; r.mode = mode; r.beta = beta;
end
