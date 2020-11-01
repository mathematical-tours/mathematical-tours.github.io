% Algorithme FFT version décimation temporelle.
function res = fft_dit(f,dir)
% N doit être une puissance de 2
N = length(f);
ldn = floor(log2(N));
f = rev_bits(f);
for ldm=1:ldn
    m = 2^ldm;
    m1 = m/2;
    for j=0:m1-1
        e = exp(-dir*2.0i*pi*j/m);
        for r=0:m:N-m
            u = f(r+j+1);
            v = f(r+j+m1+1)*e;
            f(r+j+1)    = u + v;
            f(r+j+m1+1) = u - v;
        end
    end
end
res = f; % a proscrire dans une bonne implémentation en C