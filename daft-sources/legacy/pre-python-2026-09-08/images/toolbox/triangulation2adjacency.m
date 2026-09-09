function A = triangulation2adjacency(face)

% triangulation2adjacency - compute the adjacency matrix
%   of a given triangulation.
%
%   Copyright (c) 2003 Gabriel Peyré

nvert = max(max(face));
nface = length(face);
A = zeros(nvert);

for i=1:nface
    A(face(i,1),face(i,2)) = 1;
    A(face(i,2),face(i,3)) = 1;
    A(face(i,3),face(i,1)) = 1;   
    % make sure that all edges are symmetric
    A(face(i,2),face(i,1)) = 1;
    A(face(i,3),face(i,2)) = 1;
    A(face(i,1),face(i,3)) = 1;
end