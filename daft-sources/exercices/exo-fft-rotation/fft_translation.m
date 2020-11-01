function y = fft_translation(x, v)
n = length(x);
[s,t] = meshgrid( [0:(n/2-1),-n/2:-1] );
mult = exp( 2i*pi/n*( s*v(1) + t*v(2) ) );
y = fft2(x).*mult;
y = real( ifft2(y) );