function y = resolution_implicite(x,s,theta,t)
n = length(x);  y = x;
A = zeros(n,1); A(1) = -2*s; A(2) = s; A(n) = s;
fA = fft(A); y = fft(x);
mult = ( ones(n,1)+(1-theta)*fA )./( ones(n,1)-theta*fA );
for i = 1:t
    y = y.*mult;
end
y = real( ifft(y) );