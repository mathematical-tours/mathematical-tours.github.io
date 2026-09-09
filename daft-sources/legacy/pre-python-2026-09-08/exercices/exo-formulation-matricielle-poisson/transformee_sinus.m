function y = transformee_sinus(x)
n = length(x);
x = [0;x;0;-x(n:-1:1)]; x = fft(x);
y = real( x(2:1:n+1)/(-2i) );