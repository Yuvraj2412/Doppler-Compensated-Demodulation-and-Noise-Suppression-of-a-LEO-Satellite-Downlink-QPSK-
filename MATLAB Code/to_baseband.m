function y = to_baseband(x, Fs_in, f_off, Fs_out)
% Shift signal at f_off (Hz) to 0, low-pass, resample to Fs_out (e.g. 8 x 72 kHz).
n = (0:numel(x)-1).';
x = x .* exp(-1j*2*pi*f_off*n/Fs_in);
fc = 0.42*Fs_out;  M = 255; m = (-(M-1)/2:(M-1)/2).';
h = 2*fc/Fs_in * sinc(2*fc/Fs_in*m) .* (0.54 + 0.46*cos(2*pi*m/M));
D = (M-1)/2;  x = filter(h, 1, [x; zeros(D,1)]);  x = x(D+1:end);
if abs(Fs_in - Fs_out) < 1e-6, y = x; return; end
if exist('resample','file') == 2 || exist('resample','builtin') == 5
    [P,Q] = rat(Fs_out/Fs_in, 1e-9);  y = resample(x, P, Q);
else
    t = (0:floor((numel(x)-1)*Fs_out/Fs_in)).'/Fs_out*Fs_in + 1;   % fractional indices
    y = interp1((1:numel(x)).', x, t, 'spline');
end
end
