function h = KernelDensity2D(x,n,s)

% kernel density estimator

m = 2*n; % for padding

% gaussian blur, 2D
t = [0:m/2,-m/2+1:-1]'/n;
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2))); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );

% pad = @(x)[x,zeros(m-n,1);zeros(m-n,m)];
crop = @(x)x(1:n,1:n);
GFilt = @(f,s)crop(fconv(f, G(s)));

h = zeros(n,n);
I = floor(x*(n-1))+1;

K = find(I(:,1)>0 & I(:,1)<=n & I(:,2)>0 & I(:,2)<=n);
I = I(K,:);

J = I(:,1) + (I(:,2)-1)*n;
h = hist(J,1:n*n)';
h = reshape(h, [n n]);

H = zeros(m,m); H(1:n,1:n) = h;
h = GFilt(H,s);
h = n*h; 

end


