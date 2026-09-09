% test for the dijkstra algorithm
%
%   Copyright (c) 2004 Gabriel Peyré

[M, cm] = imread('test.png');
colormap(cm);
M = double(M);

[A,vertex] = build_graph_from_image( M, 8 );
n = size(vertex,1);

gplot( A,vertex, 'k.-' );
axis tight;
axis off;


W = build_euclidean_weight_matrix(A,vertex,Inf);


start_verts = 1;  % start vertex
options.stop_at = [15+15*sqrt(n)];
options.stop_when = 100000;


% set up callback function for A*
heuristic.func  = 'heuristic_callback';
heuristic.weight  = 0;
heuristic.target = options.stop_at;
heuristic.vertex = vertex;

data = dijkstra(W, start_verts, options, heuristic);


dijkstra_plot(data, vertex, options.stop_at);
