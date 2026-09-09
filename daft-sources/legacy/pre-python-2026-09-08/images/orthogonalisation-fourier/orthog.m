function g = orthog(f)

ff = fft(f);
g = ff./abs(ff);
g = real( ifft(g) );