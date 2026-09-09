function adj_list = adjmatrix2list(A)

%   adjmatrix2list - convert from matrix adjacency representation
%       to list adjacency. The adjacency can be a weighted matrix, 
%       and linked vertices should have weight either 'Inf' or '0'.
%
%   adj_list = adjmatrix2list(A);
%
%   adj_list is a cell array of vector.
%
%   Copyright (c) 2004 Gabriel Peyré

n = size(A,1);

for i=1:n
    I = find( and( A(i,:)>0,  A(i,:)~=Inf) );
    adj_list{i} = I;
end