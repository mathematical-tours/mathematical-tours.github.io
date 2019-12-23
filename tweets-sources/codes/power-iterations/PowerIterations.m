n = 20; 
dotp = @(x,y)sum(x(:).*y(:));

E = linspace(0,1,n);
% E(end) = 1;
E(2:end)=1;
A = diag(E);

randn('seed', 123);
u = randn(n,1);

u = u/norm(u);
L = []; 
for it=1:400
v = A*u;
L(it) = dotp(u,v);
u = v/norm(v);
% L(it) = dotp(u,A*u);
end

clf;
plot(log10(abs(L-max(eig(A)))));
axis tight;

% log(e(k))= -c*k