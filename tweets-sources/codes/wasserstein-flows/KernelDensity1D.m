function h = KernelDensity1D(x,n,s)

% kernel density estimator

m = 2*n; % for padding

% gaussian blur, 2D
t = [0:m/2,-m/2+1:-1]'/n;
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2))); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );

pad = @(x)[x;zeros(m-n,1)];
crop = @(x)x(1:n);
GFilt = @(f,s)crop(fconv(pad(f), G(s)));

% gaussian blur, 1D
t = [0:m/2,-m/2+1:-1]'/n;
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2)));
fconv = @(x,y)real( ifft( fft(x).*fft(y) ) );
% GFilt = @(f,s)fconv(f, G(s));
GFilt = @(f,s)crop(fconv(pad(f), G(s)));


h = zeros(n,1);
I = floor(x*(n-1))+1;
I = I(I>0 & I<=n);
% h(I) = 1; 

h = hist(I,1:n)';

h = GFilt(h,s);
h = n*h; % /sum(h(:));

end


