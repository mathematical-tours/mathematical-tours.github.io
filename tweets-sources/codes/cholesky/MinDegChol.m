% test for minimum degree ordering

% Laplacian in 1D
n = 30;
B = [-ones(n,1),2*ones(n,1),-ones(n,1)];
S = spdiags(B,[-1 0 1],n,n);
% S(1,end)=-1; S(end,1)=-1;
% S = kron(S,eye(n)) + kron(eye(n),S);

L = chol(S,'lower');