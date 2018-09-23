function [B,beta] = ComputeBoundaries(E,F,T)

% ComputeBoundaries - compute boundary operators from simplicial homology
% 
%   [B,beta] = ComputeBoundaries(E,F,T);
%
%   beta are the betti number B(1) and B(2)
%
% You can check for boundary constraints B{k+1}*B{k}=0
%       > B{2}*B{1} % Should be 0
%       > B{3}*B{2} % Should be 0
%
%   Copyright (c) 2018 Gabriel Peyre

if nargin<3
    T = [];
end

% dimensions
n = max(E(:)); 
if isempty(E)
    n = 0;
end
N = [n, size(E,1), size(F,1), size(T,1)];
% B{1} : R^vert -> R^edge
B{1} = sparse( [1:N(2) 1:N(2)]', [E(:,1) E(:,2)], [ones(N(2),1);-ones(N(2),1)], N(2),N(1));
% B{2} : R^edge -> R^triang
I = []; J = []; V = [];
for i=1:N(3)
    f = F(i,:);
    e = {f(1:2) f(2:3) f([1 3])}; s = [1 1 -1];
    for k=1:length(e)
        I(end+1) = i;
        [~,J(end+1)] = ismember( e{k}, E, 'rows' );
        V(end+1) = s(k);
    end
end
B{2} = sparse( I,J,V, N(3),N(2));
% B{3} : R^triang -> R^tetra
I = []; J = []; V = [];
for i=1:N(4)
    f = T(i,:);
    e = {f([2 3 4]) f([1 3 4]) f([1 2 4]) f([1 2 3])}; s = [1 -1 1 -1];
    for k=1:length(e)
        I(end+1) = i;
        [~,J(end+1)] = ismember( e{k}, F, 'rows' );
        V(end+1) = s(k);
    end
end
B{3} = sparse( I,J,V, N(4),N(3));

% compute betty numbers
beta = [];
beta(1) = size(B{1},2) - rank(full(B{1})); % dim(ker(B{1})
for k=1:2
    % B{k+1} subset B{k}
    % Im(B{k}) subset Ker(B{k+1})
    % beta(k) = dim(Ker(B{k+1})) - dim(Im(B{k}))
    beta(end+1) = size(B{k+1},2) - rank(full(B{k+1})) - rank(full(B{k}));
end
beta = beta(:);

end