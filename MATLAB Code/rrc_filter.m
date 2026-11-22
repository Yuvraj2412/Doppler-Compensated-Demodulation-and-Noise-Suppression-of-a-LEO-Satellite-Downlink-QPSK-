function h = rrc_filter(beta, span, sps, tau)
% Root-raised-cosine pulse, unit energy. tau = fractional delay in samples.
if nargin < 4, tau = 0; end
t = ((-span*sps/2 : span*sps/2).' - tau) / sps;   % time in symbol periods
h = zeros(size(t));
for i = 1:numel(t)
    ti = t(i);
    if abs(ti) < 1e-9
        h(i) = 1 - beta + 4*beta/pi;
    elseif abs(abs(ti) - 1/(4*beta)) < 1e-9
        h(i) = (beta/sqrt(2))*((1+2/pi)*sin(pi/(4*beta)) + (1-2/pi)*cos(pi/(4*beta)));
    else
        h(i) = (sin(pi*ti*(1-beta)) + 4*beta*ti*cos(pi*ti*(1+beta))) / (pi*ti*(1-(4*beta*ti)^2));
    end
end
h = h / sqrt(sum(h.^2));
end
