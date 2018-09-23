function F = TriangleFromEdges(E)

m = size(E,1);
n = max(E(:));
r = size(E,2); % dimension of the simplexes

if 0
    
[K,L] = meshgrid(1:m,1:n);
U = sort([E(K(:),1),L(:)],2);
V = sort([E(K(:),2),L(:)],2);
S = ismember(U,E, 'rows') & ismember(V,E, 'rows');
I = find(S);
F = [];
for i=I(:)'
    if not(ismember(L(i), E(K(i),:)))
        F(end+1,:) = sort( [L(i), E(K(i),:)] );
    end
end

else
    
[K,L] = meshgrid(1:m,1:n);
S = ones(m*n,1); % keep it
for d=1:r
    J = 1:r; J(d) = [];
    U = sort([E(K(:),J),L(:)],2);
    S = S & ismember(U,E, 'rows');
end

I = find(S);
F = [];
for i=I(:)'
    if not(ismember(L(i), E(K(i),:)))
        F(end+1,:) = sort( [L(i), E(K(i),:)] );
    end
end


    
end

F = unique(F, 'rows');

end