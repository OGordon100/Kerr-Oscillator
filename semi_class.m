clear all;
close all;
syms x

K = -20;
dp = -20;
ep = 0;
kappa =  1;

for n = 1:1
K_value = K;
a = K_value^2;
b = 2*dp*K_value;
c = (kappa^2)/4 + dp^2;
d = -(ep^2);


p = [a b c d];

r = roots(p);
end

[~, idx] = sort(abs(imag(r)))

rnew = r(idx)