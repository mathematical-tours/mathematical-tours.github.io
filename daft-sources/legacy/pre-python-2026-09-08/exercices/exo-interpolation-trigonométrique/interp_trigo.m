function y = interp_trigo(x,eta)
N = length(x); N0 = (N-1)/2; P =N*eta;
f = fft(x);
f = eta*[f(1:N0+1); zeros(P-N,1); f(N0+2:N)];
y = real( ifft(f) );