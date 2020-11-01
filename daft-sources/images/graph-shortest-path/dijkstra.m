function data = dijkstra(W, start_verts, options, heuristic)

% dijkstra - launch the Dijkstra algorithm.
%
%   You can provide special conditions for stop in options :
%       'options.stop_at' : stop when these points are reached
%       'options.stop_when' : stop when a given number of iterations is
%          reached.
%
%   data = dijkstra(W, start_verts [, options, heuristic]);
%
%   Copyright (c) 2004 Gabriel Peyré


if nargin<3
    options.stop_at = [];
    options.stop_when = 100;
    options.verbose = 1;
end

if isfield(options,'stop_at')
    stop_at = options.stop_at;
else
    stop_at = [];
end

if isfield(options,'stop_when')
    stop_when = options.stop_when;
else
    stop_when = 100;
end

if isfield(options,'verbose')
    verbose = options.verbose;
else
    verbose = 1;
end


stop_when = min(stop_when, size(W,1));

if nargin<4
    data = dijkstra_init(W, start_verts);
else
    data = dijkstra_init(W, start_verts, heuristic);    
end

str = 'Performing Dijkstra algorithm.';
if verbose
    h = waitbar(0,str);
end

i = 0; 
while i<stop_when
    i = i+1;
    data = dijkstra_step(data);
    if verbose
        waitbar(i/stop_when, h, sprintf('Performing Dijkstra algorithm, iteration %d.', i) );
     end
    % check if we have reach one of the end points
    for j=stop_at
        if ~isempty( find(data.O==j) )
            if verbose
                close(h);
            end
            return;
        end
    end
end

if verbose
    close(h);
end