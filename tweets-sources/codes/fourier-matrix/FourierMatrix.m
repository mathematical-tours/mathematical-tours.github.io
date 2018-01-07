%%
% Vizual display of various orthobases.

rep = '../results/fourier-matrix/';
[~,~] = mkdir(rep);

addpath('../toolbox/');

upss = @(x,q,n)x(ceil(1/q:1/q:n),ceil(1/q:1/q:n),:);
ups = @(x)upss(x,round(128/size(x,1)),size(x,1));
%
resc = @(x).5 + .5*x/max(abs(x(:)));

% [J,I] = meshgrid(0:n-1,0:n-1);

% turn in color a complex matrix
Q = @(t).5 + .5*cos(2*pi*t);
T = @(X)mod(angle(X), 2*pi)/(2*pi);
A = @(X)abs(X)/max(abs(X(:)));
Rainb = @(X)cat( 3, Q(T(X)).*A(X), Q(T(X)+1/3).*A(X), Q(T(X)+2/3).*A(X) );

% real save
mysaveR = @(X,name)imwrite(resc(ups(X)), [rep name '.png'], 'png');
% complex save
mysaveC = @(X,name)imwrite(Rainb(ups(X)), [rep name '.png'], 'png');

% Fourier
F = @(n)exp(2i*pi/n*(0:n-1)'*(0:n-1));
for n=2:2:16
    mysaveC(F(n), ['fourier-' num2str(n)]);
    mysaveC(F(n)', ['ifourier-' num2str(n)]);
end

% convolution
n = 16;
r = fft(rand(n,1));
r(1)=0;
C = real( F(n)*diag(r)*F(n)' );
mysaveR(C, 'convol');
mysaveC(diag(r), 'convol-diag');

%%
% Fourier matrix factorization
n = 8;
% even/odd intersertion
I = full(sparse(1:2:n,1:n/2,ones(n/2,1),n,n) + ...
     sparse(2:2:n,n/2+1:n,ones(n/2,1),n,n));
% double Fourier
FF = zeros(n); FF(1:n/2,1:n/2) = F(n/2); FF(n/2+1:n,n/2+1:n) = F(n/2);
% diag rescale
D = diag( [ones(n/2,1); exp(2i*pi/n*(0:n/2-1)')] );
% even/odd sum
EO = [eye(n/2),eye(n/2); eye(n/2),-eye(n/2)];
% factorizationhhk
F1 = I * FF * D * EO;


mysaveR(I, 'factor-I');
mysaveR(D, 'factor-D');
mysaveR(EO, 'factor-EO');



norm(F1-F(n), 'fro')