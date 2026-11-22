function [y, fi] = suppress_interference(x, Fs, mode)
% Narrowband-interference suppression applied to the raw (pre-matched-filter) IQ.
%  'none' : bypass
%  'iir'  : 1-pole/1-zero complex notch (r = 0.998), zero on the unit circle
%  'fir'  : 511-tap linear-phase band-stop (delta minus normalised heterodyned window)
%  'nlms' : adaptive one-tap tone canceller (NLMS, reference = estimated tone)
%  'ale'  : adaptive line enhancer, 32-tap NLMS predictor, decorrelation delay
x = x(:); N = numel(x); y = x; fi = 0;
if strcmp(mode,'none'), return; end
thr = 20;                                  % detection threshold (dB): below this, no tone -> bypass
if strcmp(mode,'ale')
    fi = NaN;
else
    [fi, snr] = est_tone(x, Fs);
    if snr < thr, y = x; fi = NaN; return; end
end
w0 = 2*pi*fi/Fs; n = (0:N-1).';
switch mode
    case 'iir'
        r = 0.998;
        y = filter([1 -exp(1j*w0)], [1 -r*exp(1j*w0)], x);
    case 'fir'
        M = 511; D = (M-1)/2; m = (-D:D).';
        win = 0.54 + 0.46*cos(pi*m/D);               % Hamming
        hl = win/sum(win);                            % unit DC gain -> exact null
        h = -hl.*exp(1j*w0*m);  h(D+1) = h(D+1) + 1;
        y = filter(h, 1, [x; zeros(D,1)]);  y = y(D+1:end);   % remove group delay
    case 'nlms'
        mu = 2e-3; ref = exp(1j*w0*n); wt = 0;
        for k = 1:N
            e = x(k) - wt*ref(k);
            wt = wt + mu*e*conj(ref(k));
            y(k) = e;
        end
    case 'ale'
        M = 32; Dl = 170; mu = 0.002; wt = zeros(M,1);
        xp = [zeros(M+Dl,1); x];
        for k = 1:N
            u = xp(k+M:-1:k+1);                       % x(k-Dl) ... x(k-Dl-M+1)
            e = x(k) - wt'*u;
            wt = wt + mu*u*conj(e)/(u'*u + 1e-9);
            y(k) = e;
        end
    otherwise
        error('unknown mode');
end
end
