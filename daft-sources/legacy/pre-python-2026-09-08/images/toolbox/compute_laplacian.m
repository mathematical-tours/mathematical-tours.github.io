function lap = compute_laplacian(A)

% compute_laplacian - return the combinatorial laplacian
%   of a given adjacency matrix
%
%   lap = compute_laplacian(A);
%
%   Copyright (c) 2004 Gabriel Peyré

n = length(A);
lap = zeros(n,n);

for i=1:n
    ni = sum( A(i,:) );
    for j=1:n
        nj = sum( A(:,j) );
        % number of neighbor
        if A(i,j) ==1
            lap(i,j) = -1/sqrt(ni*nj);
        end
    end
end

for i=1:n
    lap(i,i) = -sum( lap(i,:) );
end