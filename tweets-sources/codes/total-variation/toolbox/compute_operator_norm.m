function [L,e] = compute_operator_norm(A,u, niter)

% compute_operator_norm - compute operator norm
%
%   [L,e] = compute_operator_norm(A,u,niter);
%
%   Copyright (c) 2010 Gabriel Peyre

if nargin<3
    niter = 30;
end


u = u/norm(u(:));
e = [];
for i=1:30
    v = A(u);
    e(end+1) = sum(u(:).*v(:));
    u = v/norm(v(:));
end
L = e(end);