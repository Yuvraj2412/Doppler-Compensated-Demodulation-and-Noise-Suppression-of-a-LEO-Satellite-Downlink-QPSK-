function s = qpsk_mod(bits)
% Gray-mapped QPSK, unit average symbol energy. bits: column, even length.
b1 = bits(1:2:end); b2 = bits(2:2:end);
s = ((1-2*b1) + 1j*(1-2*b2)) / sqrt(2);
end
