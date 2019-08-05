function [D,Dsvg,Ssvg] = dijkstra(W, I, options)

% perform_fm_dijstra - *slow* (matlab) implementation of Dijstra and FM
%
%   [D,Dsvg,Ssvg] = perform_fm_dijstra(W, pstart, options);
%
%   W is an (n,n) metric matrix.
%   pstart is a (2,k) starting points.
%   options.method is either 'fm' or 'dijstra'
%   
%   D is the final distance map to pstart
%   options.svg_rate gives the rate at wich Dsvg and Ssvg is filled.
%   options.niter can be used to limit the total number of steps (partial propagation). 
%   
%   Copyright (c) 2012 Gabriel Peyre

options.null = 0;
svg_rate = getoptions(options, 'svg_rate', 10);
niter = getoptions(options, 'niter', Inf);
j_tgt = getoptions(options, 'target', []);

n = size(W,1);

h = getoptions(options, 'heuristic',zeros(n,1)); 

% Initialize the distance to \(+\infty\), excepted for the boundary conditions.
D = zeros(n,1)+Inf; % current distance
D(I) = 0; 

% Initialize the state to 0 (unexplored), excepted for the boundary point to \(1\) (front).
S = zeros(n,1);
S(I) = 1; % open

iter = 0;
Dsvg = []; Ssvg = [];
while not(isempty(I)) && iter<=niter
    iter = iter+1;
    % pop from stack
    [tmp,j] = sort(D(I)+h(I)); j = j(1);
    i = I(j); I(j) = [];
    % declare dead
    S(i) = -1; 
    % The list of the neighbors
    J = find(W(i,:));
    % Remove those that are dead
    J(S(J)==-1) = [];
    % Add them to the front
    J1 = J(S(J)==0);
    I = [I(:); J1(:)];
    S(J1) = 1;
    % update neighbor values
    for j=J(:)'
        D(j) = min(D(j), D(i) + W(i,j));
    end
    % svd
    if ( svg_rate==1 || (mod(iter,svg_rate)==1) ) && nargout>2
        Dsvg(:,end+1) = D;
        Ssvg(:,end+1) = S;
    end
    if not(isempty(j_tgt)) && sum(S(j_tgt)==-1)==length(j_tgt)        
        break
    end
end

end