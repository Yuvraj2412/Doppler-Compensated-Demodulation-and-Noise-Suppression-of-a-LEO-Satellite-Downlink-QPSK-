function [u, errfrac] = viterbi_k7(r1, r2)
% Soft-decision Viterbi decoder, rate 1/2, K=7, polynomials 0x4F (first output) and 0x6D
% (second output); bit k of a polynomial taps x(n-k). Soft input: positive = coded bit 0.
% Returns decoded bits u (0/1) and the fraction of channel symbols that disagree with
% the re-encoded decoded sequence (an estimate of the raw symbol error rate).
r1 = r1(:); r2 = r2(:); N = numel(r1);
ns = (0:63).'; p0 = bitshift(ns,-1); p1 = p0 + 32; uin = bitand(ns,1);
reg0 = bitor(bitshift(p0,1), uin); reg1 = bitor(bitshift(p1,1), uin);
a0 = 1-2*par7(bitand(reg0,79)); b0 = 1-2*par7(bitand(reg0,109));
a1 = 1-2*par7(bitand(reg1,79)); b1 = 1-2*par7(bitand(reg1,109));
pm = zeros(64,1); dec = false(64,N);
for n = 1:N
    A = pm(p0+1) + a0*r1(n) + b0*r2(n);
    B = pm(p1+1) + a1*r1(n) + b1*r2(n);
    d = B > A;  pm = max(A,B);  dec(:,n) = d;  pm = pm - max(pm);
end
[~,st] = max(pm); u = zeros(N,1); pr = zeros(N,1);
for n = N:-1:1
    u(n) = bitand(st-1, 1);  p = bitshift(st-1,-1) + 32*dec(st,n);
    pr(n) = p;  st = p + 1;
end
reg = bitor(bitshift(pr,1), u);
c1 = par7(bitand(reg,79)); c2 = par7(bitand(reg,109));
errfrac = (sum((r1 < 0) ~= c1) + sum((r2 < 0) ~= c2)) / (2*N);
end

function p = par7(x)
p = zeros(size(x));
for k = 0:6, p = p + bitand(bitshift(x,-k),1); end
p = mod(p,2);
end
