function y = fft_chirp(x,g)
p = length(x); f = zeros(p-1,1); h = zeros(p-1,1);
for k=0:p-2
    j = mod(g^k,p); jj = invmod(j,p);
    f(k+1) = x(j+1); h(k+1) = exp(-2i*pi/p*jj);
end
h = ifft(fft(f).*fft(h));
y = zeros(p,1); y(1) = sum(x);
for k=0:p-2
    j = mod(g^k,p); jj = invmod(j,p);
    y(jj+1) = x(1) + h(k+1);
end