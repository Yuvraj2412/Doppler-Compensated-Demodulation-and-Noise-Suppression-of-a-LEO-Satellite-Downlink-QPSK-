function [y, fHz, th] = costas_qpsk(x, BnTs, zeta, Rs)
% 2nd-order decision-directed Costas loop for QPSK (symbol-rate input).
x = x(:);
d  = BnTs/(zeta + 1/(4*zeta));
Kp = 4*zeta*d/(1 + 2*zeta*d + d^2);
Ki = 4*d^2   /(1 + 2*zeta*d + d^2);
N = numel(x); y = zeros(N,1); fHz = zeros(N,1); th = zeros(N,1);
theta = 0; freq = 0;
for n = 1:N
    y(n) = x(n) * exp(-1j*theta);
    e = (sign(real(y(n)))*imag(y(n)) - sign(imag(y(n)))*real(y(n))) / sqrt(2);
    freq  = freq + Ki*e;
    theta = theta + freq + Kp*e;
    fHz(n) = freq/(2*pi)*Rs;  th(n) = theta;
end
end
