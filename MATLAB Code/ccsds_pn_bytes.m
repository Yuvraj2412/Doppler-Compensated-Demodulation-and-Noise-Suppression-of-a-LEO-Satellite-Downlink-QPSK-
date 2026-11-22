function b = ccsds_pn_bytes(nBytes)
% CCSDS pseudo-randomiser, h(x) = x^8+x^7+x^5+x^3+1, all-ones seed.
% Recurrence a(n+8) = a(n+7) ^ a(n+5) ^ a(n+3) ^ a(n). Starts FF 48 0E C0 9A 0D 70 BC.
nb = 8*nBytes; a = zeros(1, nb+8); a(1:8) = 1;
for n = 1:nb, a(n+8) = xor(xor(a(n+7), a(n+5)), xor(a(n+3), a(n))); end
b = (reshape(a(1:nb), 8, []).' * 2.^(7:-1:0).').';
end
