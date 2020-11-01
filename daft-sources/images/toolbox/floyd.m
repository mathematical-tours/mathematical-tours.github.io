function D1 = floyd(D, verbose)

% floyd - compute shortest distance on a 
%   graph using Floyd algorithm.
%   The matrix A is a weighted adjacency matrix
%   The weight should be Inf if the 2 vertices
%   are not connected.
%
%   D1 = floyd(D, verbose);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin==1
    verbose = 1;
end

N = length(D);

if verbose
    h = waitbar(0,'Computing shortest path distance.');
end

for k=1:N
     D = min(D,repmat(D(:,k),[1 N])+repmat(D(k,:),[N 1])); 
     if verbose
         waitbar(k/N);
     end
end

if verbose
    close(h);
end

D1 = D;