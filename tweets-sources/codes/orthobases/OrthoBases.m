%%
% Vizual display of various orthobases.

rep = '../results/orthobases/';
[~,~] = mkdir(rep);

addpath('../toolbox/');

n = 16;

q= round(128/n);
ups = @(x)x(ceil(1/q:1/q:n),ceil(1/q:1/q:n));
resc = @(x).5 + .5*x/max(abs(x(:)));

[J,I] = meshgrid(0:n-1,0:n-1);

mysave = @(X,name)imwrite(resc(ups(X)), [rep name '.png']);

% identity
X = eye(n);
mysave(X, 'identity');

% walsh
X0 = [1 1;-1 1]/sqrt(2);
X = X0;
for j=1:log2(n)-1
    X = kron(X,X0);
end
norm(X*X'-eye(n), 'fro')
mysave(X,  'walsh');

% dct2
X = cos(pi/n*(I+1/2).*J); X(:,1) = X(:,1)*1/sqrt(2); X = X*sqrt(2/n);
norm(X*X'-eye(n), 'fro')
mysave(X,  'dct2');

% dct4
X = cos(pi/n*(I+1/2).*(J+1/2)); X = X*sqrt(2/n);
norm(X*X'-eye(n), 'fro')
mysave(X,  'dct4');


% random orthogonal
[X,R] = qr(randn(n));
norm(X*X'-eye(n), 'fro')
mysave(X, 'rand');


% random permutation
X = eye(n); X = X(randperm(n),:);
norm(X*X'-eye(n), 'fro')
mysave(X, 'randperm');

% haar
for i=1:n
    a = zeros(n,1); a(i)=1;
    X(:,i) = perform_haar_transf(a,0,-1);
end
norm(X*X'-eye(n), 'fro')
mysave(X,  'haar');

% random orthogonal convolution
U = exp(2i*pi/n*I.*J);
u = fft(randn(n,1)); u = u./abs(u);
X = real( U*diag(u)*U'/n );
norm(X*X'-eye(n), 'fro')
mysave(X,  'randconv');

% daub4
options.h =  compute_wavelet_filter('Daubechies',4);
for i=1:n
    a = zeros(n,1); a(i)=1;
    X(:,i) = perform_wavortho_transf(a,1,-1,options);
end
norm(X*X'-eye(n), 'fro')
mysave(X, 'daub4');