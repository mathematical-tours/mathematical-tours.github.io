function W = perform_cart(T)

% perform_cart - perform CART algorithm
%
%   W = perform_cart(T);
%
%   T is a weight tree, where T{j}(i) is the cost on node i at depth j.
%   W is a tree where W{j}(i)=0 for an interior node, +1 for a leaf, and -1
%       for a non used node (below a leaf).
%
%   T is assumed to be a uniform tree with a fixed branching factor.
%
%   W is the tree that minimize sum_{(i,j) leaf of W} T{j}(i).
%
%   Copyright (c) 2010 Gabriel Peyre

% branching factor
q = length(T{2}/T{1});
% depth
J = length(T);


%%
% Step 1: build the tree
d = T{J};
W = {}; W{J} = ones(q^(J-1),1); % leaves
for j=J-1:-1:1
    % leaf contributions
    d1 = sum(reshape(d, [q length(d)/q]))';
    I = find(d1<=T{j});
    d = min(d1, T{j});
    W{j} = ones(q^(j-1),1);
    W{j}(d~=T{j}) = 0; % interior nodes
%    W{j}(d==T{j}) = 1; % new leaves
end
%%
% Step 2: prune the tree (set to -1 dead nodes)
for j=2:J
    I = find(W{j-1}==-1 | W{j-1}==1);
    for s=1:q
        W{j}(q*(I-1)+s) = -1;
    end
end
