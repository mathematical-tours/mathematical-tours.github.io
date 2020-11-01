% Algorithme FFT version décimation fréquentielle.
function res = fft_dif(f,dir)
% N doit être une puissance de 2
N = length(f);
ldn = floor(log2(N));
for ldm=ldn:-1:1
    m = 2^ldm;
    m1 = m/2;
    for j=0:m1-1
        e = exp(-dir*2.0i*pi*j/m);
        for r=0:m:N-1
            u = f(r+j+1);
            v = f(r+j+m1+1);
            f(r+j+1)    = u + v;
            f(r+j+m1+1) = (u-v)*e;
        end
    end
end
% remet le vecteur transformé dans le bon ordre
f = rev_bits(f);
res = f;  % a proscrire dans une bonne implémentation en C