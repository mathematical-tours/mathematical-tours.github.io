function y = fft_transvec_x(x, lambda)
n = length(x);
for k=1:n
    v = x(:,k); trans = lambda*(k-n/2-1);
    mult = exp( 2i*pi/n*( [0:(n/2-1),-n/2:-1]'*trans ) );
    v = fft(v).*mult; y(:,k) = real( ifft(v) );
end