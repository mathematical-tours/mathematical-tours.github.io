function y = resolution_implicite(x,s,theta,t)
n = length(x);  y = x;
A = zeros(n,n); A(1,1)=-4*s; A(2,1)=s; A(1,2)=s; A(n,1)=s; A(1,n)=s;
fA = fft2(A); y = fft2(x);
mult = ( ones(n,n)+(1-theta)*fA )./( ones(n,n)-theta*fA );
for i = 1:t
    y = y.*mult;
end
y = real( ifft2(y) );