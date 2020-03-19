function [L,W,D,Di] = MeshLaplacian(X,F, lapl)

if nargin<3
    lapl = 'cotan';
end

n = size(X,2);
W = sparse(n,n);

switch lapl
    case 'cotan'
        for i=1:3
            i2 = mod(i  ,3)+1;
            i3 = mod(i+1,3)+1;
            u = X(:,F(i2,:)) - X(:,F(i,:));
            v = X(:,F(i3,:)) - X(:,F(i,:));
            % normalize the vectors
            u = u ./ repmat( sqrt(sum(u.^2,1)), [3 1] );
            v = v ./ repmat( sqrt(sum(v.^2,1)), [3 1] );
            % compute angles
            alpha = 1 ./ tan( acos(sum(u.*v,1)) );
            alpha = max(alpha, 1e-2); % avoid degeneracy
            W = W + sparse(F(i2,:),F(i3,:), alpha, n, n );
            W = W + sparse(F(i3,:),F(i2,:), alpha, n, n );
        end
    case 'combinatorial'
        E = [F([1 2],:) F([2 3],:) F([3 1],:)];
        p = size(E,2);
        W = sparse( E(1,:), E(2,:), ones(p,1) );
        W = max(W,W');
        
end

d = full( sum(W,1) );
D = spdiags(d(:), 0, n,n);
Di = spdiags(1./d(:), 0, n,n);
L = D - W;

end