function out = run_link_oqpsk(EbN0_dB, nSym, cfg)
% OQPSK version of run_link: TX (Q delayed T/2) -> Doppler S-curve + clock offset
%  (+tone) -> AWGN -> [suppression] -> RRC MF -> coarse CFO (4th power)
%  -> joint OQPSK timing/carrier loop -> BER (I and Q branches scored separately).
sps = cfg.sps; Rs = cfg.Rs; Fs = Rs*sps; L = cfg.span*sps + 1;
hT = rrc_filter(cfg.beta, cfg.span, sps, cfg.tauSym*sps);
h  = rrc_filter(cfg.beta, cfg.span, sps);

bits = randi([0 1], 2*nSym, 1);  s = qpsk_mod(bits);
up = zeros(nSym*sps,1);
up(1:sps:end) = real(s);                          % I on the symbol grid
up(1+sps/2:sps:end) = up(1+sps/2:sps:end) + 1j*imag(s);   % Q delayed by T/2
tx = filter(hT, 1, [up; zeros(L-1,1)]);

n = (0:numel(tx)-1).';  t = n/Fs;  T = t(end);
fd = -cfg.fmax*tanh(4*(t/T - 0.5));
rx = tx .* exp(1j*2*pi*cumsum(fd)/Fs);
if cfg.ppm ~= 0, rx = interp1(n, rx, n*(1 + cfg.ppm*1e-6), 'spline', 0); end
if isfield(cfg,'intf') && ~isempty(cfg.intf)
    Ps = 1/sps;  A = sqrt(Ps/10^(cfg.intf.sir_dB/10));
    rx = rx + A*exp(1j*(2*pi*cfg.intf.f*n/Fs + 2*pi*rand));
end
sigma2 = 0.5/10^(EbN0_dB/10);
rx = rx + sqrt(sigma2/2)*(randn(size(rx)) + 1j*randn(size(rx)));

mode = 'none'; if isfield(cfg,'supp'), mode = cfg.supp; end
[rx, out.fi] = suppress_interference(rx, Fs, mode);

m = filter(h, 1, [rx; zeros(L-1,1)]);
Nf = 16384; NF = 2^18;  seg = m(8193:8192+Nf);
X = abs(fft(seg.^4, NF)); [~,k] = max(X);
f4 = (k-1)/NF*Fs; if f4 > Fs/2, f4 = f4 - Fs; end
f0 = f4/4;
m = m .* exp(-1j*2*pi*f0*(0:numel(m)-1).'/Fs);

bw = [0.005 0.01]; if isfield(cfg,'loopBW'), bw = cfg.loopBW; end   % [timing carrier] BnT
[yc, fHz] = oqpsk_sync(m, sps, bw(1), bw(2), 0.707, Rs);

% BER: 4 rotations x independent lag for I and Q branches (data-aided, sim only)
skip = 3000; b1 = bits(1:2:end); b2 = bits(2:2:end); best = inf;
idx = skip+1 : min(numel(yc), nSym-skip);
for r = 0:3
    z = yc(idx) * 1j^r;  eI = inf; eQ = inf;
    for d = -12:12
        eI = min(eI, sum((real(z) < 0) ~= b1(idx+d)));
        eQ = min(eQ, sum((imag(z) < 0) ~= b2(idx+d)));
    end
    best = min(best, eI + eQ);
end
out.nbits = 2*numel(idx); out.ber = best/out.nbits;
out.y = yc; out.fest = f0 + fHz; out.ftrue = fd(1:sps:end); out.f0 = f0;
end
