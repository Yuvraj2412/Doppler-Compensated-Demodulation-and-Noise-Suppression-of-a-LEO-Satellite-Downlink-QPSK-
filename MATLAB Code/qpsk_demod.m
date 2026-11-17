function bits = qpsk_demod(s)
s = s(:);
bits = zeros(2*numel(s),1);
bits(1:2:end) = real(s) < 0;
bits(2:2:end) = imag(s) < 0;
end
