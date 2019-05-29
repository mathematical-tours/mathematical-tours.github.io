%%
% Display weighted laplacian eigenvectors.

n = 30; 


% grad operator with neumann BC
dx = spdiags([-ones(n,1), ones(n,1)], -1:0, n, n);
dx = dx(2:end,:);

% with periodic BC
dx = spdiags([-ones(n,1), ones(n,1)], -1:0, n, n);
% dx(1,end) = -1;


n1 = size(dx,1);

t = linspace(0,1,n1)';
[Y,X] = meshgrid(t);
gauss = @(m,s)exp(-(X-m(1)).^2/(2*s).^2-(Y-m(2)).^2/(2*s).^2);

Grad = [kron(dx,speye(n)); kron(speye(n),dx)];


% two weight functions
a = ones(n1);

b = gauss([.2,.25],.05) + gauss([.7,.74],.06);
% b = 1-b;
b = b + .02;

K = 50; %# eigenvect 

q = 50;
V1 = eye(size(Grad,2),K);
for it=1:q
    s = (it-1)/(q-1);
    c = (1-s)*a+s*b; 
    Delta = Grad' * spdiags( [c(:); c(:)], 0, 2*n*n,2*n*n ) * Grad;   
    [V,D] = eigs(Delta,K, 'SM');
    V = V * diag(diag(sign(V'*V1)));
    V1 = V;
    %
    v = reshape(V(:,40), [n n]);
    clf; imagesc(v);
    axis tight;
    drawnow;
end