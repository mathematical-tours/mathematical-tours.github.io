function data = dijkstra_init(W, start_verts, heuristic)

% dijkstra_init - initialisation of dijkstra algorithm
%
%   data = dijkstra_init(W, start_verts [,heuristic]);
%
%   'heuristic' is a structure that should contains :
%       - one field 'heuristic.func' which should be a function.
%           This function take as input the number of a vertex
%           and should return a heuristical measure of the distance
%           between this point and the target.
%       - one field 'heuristic.weight' in [0,1] which measure
%           how much this heuristic should be taken into acount.
%           0 is classical Dijkstra, and 1 is full A* algorithm.
%       - you can add other fields (use data) for your function.
%
%   Copyright (c) 2004 Gabriel Peyré

n = size(W,1);

if nargin<3
    heuristic.func  = 0;
    heuristic.weight  = 0;
end

data.heuristic = heuristic;

data.A = zeros(n,1) + Inf; % action 
data.A(start_verts) = 0;

data.O = start_verts;

data.C = [];

data.F = zeros(n,1) - 1;
data.H = zeros(n,1);
data.S = zeros(n,1);
data.S(start_verts) = 'O';

data.adj_list = adjmatrix2list(W);

data.W = W;