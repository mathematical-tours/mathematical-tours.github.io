% just quadratic functions

% quad matrix
t = .3; 
mu = .2;
U = [cos(t) sin(t);-sin(t) cos(t)];
A = U*diag([1,mu])*U';
Ai = inv(A);

% display grid
q = 256;
t = linspace(-1,1,q);
[Y,X] = meshgrid(t,t);
F = A(1,1)*X.^2+A(2,2)*Y.^2+2*A(1,2).*X.*Y;