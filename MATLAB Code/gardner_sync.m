function y = gardner_sync(x, sps, BnTs, zeta)
% Gardner timing recovery with cubic (Lagrange) interpolation.
% x: matched-filter output at sps samples/symbol. Returns 1 sample/symbol.
x = x(:);  Kted = 2.7;
d  = BnTs/(zeta + 1/(4*zeta));
Kp = 4*zeta*d/(1 + 2*zeta*d + d^2)/Kted;
Ki = 4*d^2   /(1 + 2*zeta*d + d^2)/Kted;
ip = @(p) interp4(x, p);
pos = sps + 2;  step = sps;  yprev = ip(pos);  intg = 0;
y = zeros(floor(numel(x)/sps),1); n = 0;
while pos + step + 3 < numel(x)
    pos  = pos + step;
    ycur = ip(pos);  ymid = ip(pos - step/2);
    e = real((ycur - yprev)*conj(ymid));
    e = max(min(e, 1), -1);
    intg = intg + Ki*e;
    step = sps*(1 - (Kp*e + intg));     % late (e>0) -> shorten interval
    n = n + 1; y(n) = ycur; yprev = ycur;
end
y = y(1:n);
end

function v = interp4(x, p)
i = floor(p); f = p - i;
c = [-f*(f-1)*(f-2)/6, (f+1)*(f-1)*(f-2)/2, -(f+1)*f*(f-2)/2, (f+1)*f*(f-1)/6];
v = c * x(i-1:i+2);
end
