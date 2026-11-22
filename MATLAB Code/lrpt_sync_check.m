function R = lrpt_sync_check(y)
% Tries the 8 QPSK symbol-mapping variants (4 rotations x I/Q swap), Viterbi-decodes
% each, then looks for the 32-bit sync marker 0x1ACFFC1D (or its complement) and checks
% that it repeats every 8192 bits (1024-byte frame). Also reports the re-encode
% symbol-error fraction: ~11-13 % for a wrong variant, near 0 for the right one (a 180-degree
% rotation also decodes cleanly, to the complemented bit stream, hence the complement search).
y = y(:);
asm = double(dec2bin(hex2dec('1ACFFC1D'), 32).' == '1');
R = struct('sw',{},'rot',{},'err',{},'best',{},'inv',{},'hits',{},'periodic',{},'M',{},'u',{});
v = 0;
for sw = 0:1
    for rot = 0:3
        z = y * (1j)^rot;  if sw, z = 1j*conj(z); end      % swap I and Q
        [u, err] = viterbi_k7(real(z), imag(z));
        M = zeros(numel(u)-31, 1);
        for k = 1:32, M = M + (u(k:end-32+k) == asm(k)); end
        Mc = 32 - M;  Mb = max(M, Mc);
        v = v + 1;  R(v).sw = sw; R(v).rot = rot; R(v).err = err; R(v).M = Mb; R(v).u = uint8(u);
        R(v).best = max(Mb); R(v).inv = max(Mc) > max(M);
        hit = find(Mb >= 30);  R(v).hits = numel(hit);     % <=2 bit errors: chance hits ~0.02 per variant
        per = 0; for h = hit.', per = per + any(abs(hit - (h + 8192)) <= 1); end
        R(v).periodic = per;
    end
end
end
