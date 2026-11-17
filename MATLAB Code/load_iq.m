function [x, Fs] = load_iq(file, fmt, Fs, skipSec, durSec)
% Read raw IQ. fmt: 'cu8' (RTL-SDR), 'cs8', 'cs16', 'cf32' (interleaved I,Q) or 'wav' (2-ch).
if strcmp(fmt,'wav')
    info = audioinfo(file);  Fs = info.SampleRate;
    a = round(skipSec*Fs) + 1;  b = min(info.TotalSamples, a + round(durSec*Fs) - 1);
    y = audioread(file, [a b]);                    % read only the requested segment
    x = complex(y(:,1), y(:,2)); return;
end
switch fmt
    case 'cu8',  typ = 'uint8';   nb = 1;
    case 'cs8',  typ = 'int8';    nb = 1;
    case 'cs16', typ = 'int16';   nb = 2;
    case 'cf32', typ = 'float32'; nb = 4;
    otherwise, error('unknown format');
end
fid = fopen(file,'rb'); if fid < 0, error('cannot open %s', file); end
fseek(fid, round(skipSec*Fs)*2*nb, 'bof');
raw = double(fread(fid, 2*round(durSec*Fs), [typ '=>' typ])); fclose(fid);
if strcmp(fmt,'cu8'), raw = raw - 127.5; end
raw = raw(1:2*floor(numel(raw)/2));
x = complex(raw(1:2:end), raw(2:2:end));
end
