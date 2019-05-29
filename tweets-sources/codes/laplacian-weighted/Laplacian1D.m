%%
% Display weighted laplacian eigenvectors.

n = 256*4; 
t = linspace(0,1,n)';
gauss = @(m,s)exp(-(t-m).^2/(2*s).^2);


% periodic FD
dx = spdiags([-ones(n,1), ones(n,1)], -1:0, n, n);

dx = dx(2:end,:);
% dx(1,end) = -1;

dx = full(dx);
% two weight functions
a = ones(n,1);


s = .6;
b = .05 + gauss(.3*s,.05*s) + .5*gauss(.6*s,.15*s);
b = [ones(n/2,1);ones(n/2,1)*30];
% b = .5*gauss(1-s*.2-t,.04*s) + .8*gauss(1-s*.5-t,.05*s) + .5*gauss(1-s*.8-t,.04*s);
b = gauss(.5,.08) + .05;

b = gauss(.2,.03) + .9*gauss(.7,.05) + .02;

q = 50;
V1 = eye(size(dx,1));
for it=1:q
    s = (it-1)/(q-1);
    c = (1-s)*a+s*b;
    M = dx*diag(1./c)*dx';
    [V,D] = eig(M);
    V = V * diag(diag(sign(V'*V1)));
    V1 = V;
    clf; plot(V(:,50));
    axis tight;
    drawnow;
end