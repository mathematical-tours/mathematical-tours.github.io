function d = heuristic_callback(i, heuristic)

% heuristic_callback - a default callback for A*
%
%   d = heuristic_callback(i, heuristic);
%
%   'heuristic' should contains a field 
%   'target' a target vertex number and
%   'vertex' a matrix with position of all the vertices.
%
%   Copyright (c) 2004 Gabriel Peyré

vertex = heuristic.vertex;
target = heuristic.target;

d = norm( vertex(i,:)-vertex(target,:)  , 'fro' );