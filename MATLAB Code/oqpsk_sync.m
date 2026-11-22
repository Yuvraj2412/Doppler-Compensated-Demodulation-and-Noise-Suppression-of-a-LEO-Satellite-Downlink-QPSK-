function [y, fHz] = oqpsk_sync(x, sps, BnT_t, BnT_c, zeta, Rs)
% Joint timing + carrier recovery for OQPSK (symbol-synchronous loop).
% x: matched-filter output, sps samples/symbol (sps even), residual CFO small.
% At each symbol k (strobe pos p): I peaks at p, Q peaks at p + T/2.
%   timing error : Gardner per PAM branch (I: mid=p-T/2, Q: mid=p)
%   carrier error: decision-directed QPSK Costas on (I_k, Q_k)
% Returns y(k) = I_k + j*Q_k (QPSK-like, after rotation by the carrier loop).
x = x(:);  Kted = 2.7;
dt = BnT_t/(zeta + 1/(4*zeta));
Kpt = 4*zeta*dt/(1+2*zeta*dt+dt^2)/Kted;  Kit = 4*dt^2/(1+2*zeta*dt+dt^2)/Kted;
dc = BnT_c/(zeta + 1/(4*zeta));
Kpc = 4*zeta*dc/(1+2*zeta*dc+dc^2);       Kic = 4*dc^2/(1+2*zeta*dc+dc^2);
ip = @(p) interp4(x, p);
pos = 2*sps + 2;  step = sps;  theta = 0;  freq = 0;  intg = 0;
Iprev = 0; Qprev = 0;
N = floor(numel(x)/sps); y = zeros(N,1); fHz = zeros(N,1); n = 0;
while pos + 1.6*step + 4 < numel(x)
    pos = pos + step;
    rot = exp(-1j*theta);
    a = ip(pos - step/2)*rot;  b = ip(pos)*rot;  c = ip(pos + step/2)*rot;
    Ik = real(b);  Qk = imag(c);
    e  = real(a)*(Ik - Iprev) + imag(b)*(Qk - Qprev);     % late -> e > 0
    e  = max(min(e, 1), -1);
    intg = intg + Kit*e;
    step = sps*(1 - (Kpt*e + intg));
    ec = (sign(Ik)*Qk - sign(Qk)*Ik)/sqrt(2);
    freq  = freq + Kic*ec;
    theta = theta + freq + Kpc*ec;
    Iprev = Ik;  Qprev = Qk;
    n = n + 1;  y(n) = Ik + 1j*Qk;  fHz(n) = freq/(2*pi)*Rs;
end
y = y(1:n); fHz = fHz(1:n);
end

function v = interp4(x, p)
i = floor(p); f = p - i;
c = [-f*(f-1)*(f-2)/6, (f+1)*(f-1)*(f-2)/2, -(f+1)*f*(f-2)/2, (f+1)*f*(f-1)/6];
v = c * x(i-1:i+2);
end
