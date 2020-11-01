function y = dct3(x)
n = length(x);
y = [x ; zeros(3*n,1)];
y = real( fft(y) );
y = y(2:2:2*n) - x(1)/2;