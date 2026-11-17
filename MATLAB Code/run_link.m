function out = run_link(EbN0_dB, nSym, cfg)
% Simulated chain: TX -> Doppler S-curve + clock offset -> (+interferer) -> AWGN
%   -> interference suppression -> RRC MF -> coarse CFO -> Gardner -> Costas -> BER
% Optional cfg fields: cfg.intf.f (Hz), cfg.intf.sir_dB ; cfg.supp (mode string)
sps = cfg.sps; Rs = cfg.Rs; Fs = Rs*sps; L = cfg.span*sps + 1;
hT = rrc_filter(cfg.beta, cfg.span, sps, cfg.tauSym*sps);
h  = rrc_filter(cfg.beta, cfg.span, sps);

bits = randi([0 1], 2*nSym, 1);  s = qpsk_mod(bits);
up = zeros(nSym*sps,1);  up(1:sps:end) = s;
tx = filter(hT, 1, [up; zeros(L-1,1)]);

n = (0:numel(tx)-1).';  t = n/Fs;  T = t(end);
fd = -cfg.fmax*tanh(4*(t/T - 0.5));
rx = tx .* exp(1j*2*pi*cumsum(fd)/Fs);
if cfg.ppm ~= 0
    rx = interp1(n, rx, n*(1 + cfg.ppm*1e-6), 'spline', 0);
end
if isfield(cfg,'intf') && ~isempty(cfg.intf)          % fixed-frequency ground interferer
    Ps = 1/sps;  A = sqrt(Ps/10^(cfg.intf.sir_dB/10));
    rx = rx + A*exp(1j*(2*pi*cfg.intf.f*n/Fs + 2*pi*rand));
end
sigma2 = 0.5/10^(EbN0_dB/10);
rx = rx + sqrt(sigma2/2)*(randn(size(rx)) + 1j*randn(size(rx)));
out.rx_pre = rx(1:min(numel(rx),2^17));

mode = 'none'; if isfield(cfg,'supp'), mode = cfg.supp; end
[rx, out.fi] = suppress_interference(rx, Fs, mode);
out.rx_post = rx(1:min(numel(rx),2^17));

m = filter(h, 1, [rx; zeros(L-1,1)]);
Nf = 16384; NF = 2^18;
seg = m(8193:8192+Nf);   % skip filter-convergence transients
X = abs(fft(seg.^4, NF)); [~,k] = max(X);
f4 = (k-1)/NF*Fs; if f4 > Fs/2, f4 = f4 - Fs; end
f0 = f4/4;
m = m .* exp(-1j*2*pi*f0*(0:numel(m)-1).'/Fs);

ys = gardner_sync(m, sps, 0.005, 0.707);
[yc, fHz] = costas_qpsk(ys, 0.01, 0.707, Rs);

skip = 3000; best = 1; nb = 1;
for r = 0:3
  for d = -12:12
    idx = skip+1 : min(numel(yc), nSym-skip);
    ref = s(idx + d);  rec = yc(idx) * 1j^r;
    errs = sum(qpsk_demod(rec) ~= qpsk_demod(ref));
    if errs/(2*numel(idx)) < best, best = errs/(2*numel(idx)); nb = 2*numel(idx); end
  end
end
out.ber = best; out.nbits = nb; out.y_pre = ys; out.y = yc;
out.fest = f0 + fHz; out.ftrue = fd(1:sps:end); out.f0 = f0;
end
