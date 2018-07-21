function h = KernelDensity1D(x,n,s)

% kernel density estimator

% gaussian blur, 2D
t = [0:n/2,-n/2+1:-1]'/n;
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2))); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));

% gaussian blur, 1D
t = [0:n/2,-n/2+1:-1]'/n;
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2)));
fconv = @(x,y)real( ifft( fft(x).*fft(y) ) );
GFilt = @(f,s)fconv(f, G(s));


h = zeros(n,1);
I = floor(x*(n-1))+1;
% h(I) = 1; 

h = hist(I,1:n)';

h = GFilt(h,s);
h = n*h/sum(h(:));

end


