% STAGE 4a: load a real LEO IQ recording, find the signal, bring it to a
% 576 kHz (= 8 x 72 ksym/s) complex baseband file, and clean up tone interferers.
% EDIT THE CONFIG BLOCK, run, check the plots, then use the saved .mat in Stage 4b.
clear; clc; close all;
%% ---- CONFIG -------------------------------------------------------------
file = 'C:\Users\tomar\Downloads\DSP_Project\LRPT_16-29-02_137900kHz.wav';
fmt = 'wav';   Rs = 72e3;   skipSec = 0;   durSec = 10;
Fs_in   = 1.024e6;                % sample rate of the file (Hz); ignored for wav
skipSec = 0;  durSec = 10;        % portion to process
f_off   = [];                     % signal centre offset in file (Hz); [] = auto-estimate
Rs = 72e3; sps = 8; Fs_out = Rs*sps;
%% -------------------------------------------------------------------------
[x, Fs_in] = load_iq(file, fmt, Fs_in, skipSec, durSec);
x = x - mean(x);                                       % remove DC spike (RTL-SDR)
fprintf('Loaded %d samples (%.1f s) at %.0f Sa/s\n', numel(x), numel(x)/Fs_in, Fs_in);

% Overview spectrum / spectrogram of the raw file
nseg = min(numel(x), round(5*Fs_in));
[S, f, t] = simple_spectrogram(x(1:nseg), Fs_in, 4096, 4096);
P = mean(10.^(S/10), 2);
figure; subplot(2,1,1); plot(f/1e3, 10*log10(P)); grid on; xlabel('kHz'); ylabel('dB/Hz'); title('Raw PSD');
subplot(2,1,2); imagesc(t, f/1e3, S); axis xy; colorbar; xlabel('Time (s)'); ylabel('kHz'); title('Raw spectrogram');

% Auto-estimate centre/bandwidth of the wide signal (half-power edges above the floor)
Pp = [P(end-30:end); P; P(1:31)]; Ps = zeros(size(P));      % running median (rejects narrow spikes/tones)
for i = 1:numel(P), Ps(i) = median(Pp(i:i+60)); end
Ps = filter(ones(21,1)/21, 1, Ps); Ps = circshift(Ps, -10);
Pss = sort(Ps);  fl = median(Pss(1:round(0.2*numel(Pss))));  pk = max(Ps);
above = find(Ps > fl + 0.5*(pk - fl));
f_lo = f(above(1)); f_hi = f(above(end));
fprintf('Estimated signal band: %.1f .. %.1f kHz  (centre %.1f kHz, width %.1f kHz)\n', ...
        f_lo/1e3, f_hi/1e3, (f_lo+f_hi)/2e3, (f_hi-f_lo)/1e3);
if isempty(f_off), f_off = (f_lo+f_hi)/2; end
fprintf('Using f_off = %.1f Hz\n', f_off);

% To baseband @ 576 kHz
xb = to_baseband(x, Fs_in, f_off, Fs_out);
xb = xb / sqrt(mean(abs(xb).^2));                      % unit RMS
[S2, f2, t2] = simple_spectrogram(xb, Fs_out, 1024, 512);
figure; subplot(1,2,1); plot(f2/1e3, 10*log10(mean(10.^(S2/10),2))); grid on
xlabel('kHz'); ylabel('dB'); title('Baseband PSD (576 kS/s)');
subplot(1,2,2); imagesc(t2, f2/1e3, S2); axis xy; colorbar; xlabel('Time (s)'); ylabel('kHz');
title('Baseband spectrogram (watch the Doppler drift)');

% Interference check (same detector as Stage 3)
[fi, snr] = est_tone(xb, Fs_out);
fprintf('Strongest narrowband line: %.1f Hz, %.1f dB above local level -> ', fi, snr);
if snr >= 20, fprintf('interferer, applying IIR notch\n'); xb = suppress_interference(xb, Fs_out, 'iir');
else, fprintf('none detected (bypass)\n'); end

save('meteor_bb.mat', 'xb', 'Fs_out', 'Rs', 'sps', 'f_off');
fprintf('Saved meteor_bb.mat for Stage 4b.\n');
