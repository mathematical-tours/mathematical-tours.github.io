function [Q,tx,ty] = Julia(c,z,a,n,niter)

% c: fractal parameter z^2+c
% z: center
% a: radius



t = linspace(-a,a,n);
tx = real(z)+t;
ty = imag(z)+t;
[Y,X] = meshgrid(ty,tx);
Z0 = X+1i*Y;

Z = Z0;
Q = zeros(n)+Inf;
rmax = 2;
for i=1:niter
    Z = Z.^2+c;
    Q(abs(Z)>=rmax) = min(Q(abs(Z)>=rmax), i);
end
Q(Q==Inf)=niter;

end