function y = dct2(x)
n = length(x);
y(2:2:2*n,1) = x;
y = [y ; zeros(2*n,1)];
y = real( fft(y) );
y = y(1:n);