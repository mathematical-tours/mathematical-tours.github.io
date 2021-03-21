%%
% Quantum oscillator eigenfunctions

n = 256;
t = 2*(0:n-1)'/n-1;

% F'=F^{-1}
[Y,X] = meshgrid(0:n-1,0:n-1);
F = 1/sqrt(n) * exp(2i*pi/n * X.*Y);

% d^2/dx^2
D2x = ( (mod(X-Y-1,n)==0) + (mod(X-Y+1,n)==0) - 2*(X==Y) );
Q = D2x - 1000*diag(t.^2);

[U,D] = eig(Q);
plot(U(:,end-3:end));


    