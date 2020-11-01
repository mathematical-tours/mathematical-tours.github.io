function y = der_frac(f, alpha)
n = length(f)/2; 
fce = pi*[(0:n)/n, ((-n+1):-1)/n];
f = fft( f );
f = (-i*fce).^alpha.*f;
f = real(ifft(f))