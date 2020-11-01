% multiplication par convolution acyclique
function r = mult_convol(p,q)
n = length(p);
p = [p;zeros(n-1,1)];
q = [q;zeros(n-1,1)];
r = real( ifft(fft(p).*fft(q)) );