addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it,name)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


n = 2000;
n = 1000;

name = 'helix';
name = 'swiss';
noise = 0.1;
[X, labels, t] = generate_data(name, n, noise);
X(:,2) = X(:,2)*2;



clf;
scatter3(X(:,1), X(:,2), X(:,3), 15, t(:,1), 'filled');
axis equal;
axis off;
view(-20,20);

k = 5;

d = sum(X.^2,2);
D = max(d + d' - 2*X*X',0);

[Ds,Is] = sort(D, 2);
I = Is(:,2:k+1);
J = repmat((1:n)', [1 k]);
V = ones(n,k);
K = sparse(J(:),I(:),V(:), n,n);
K = max(K,K');

if 0
    sigma = mean(sqrt(D(:)))*.1;
    K = exp(-D/(2*sigma^2));
    K = double( D<=sigma^2 );
end



if 0
    L = eye(n) - K./sum(K,2);
    [U,Delta] = eig(L); Delta = real(diag(Delta));
    [Delta,I] = sort(Delta, 'ascend'); U = real(U(:,I));
else
    L = speye(n) - K./sum(K,2);
    [U,Delta] = eigs(L,5,'smallestabs');
end



U(:,2) = U(:,2) * sign(sum(U(:,2).*t(:,1)));
U(:,3) = U(:,3) * sign(sum(U(:,3).*t(:,2)));
clf;
scatter(U(:,2), U(:,3), 25, t(:,1), 'filled');
axis equal; axis off;
mysaveas(it,'flatten');

clf;
scatter3(X(:,1), X(:,2), X(:,3), 25, U(:,2), 'filled');
view(-20,20)
axis equal; axis off; axis vis3d;
q = 100;
for it=1:q
    view((it-1)/q*360,20); drawnow;
    % mysaveas(it,'eig1');
end


clf;
scatter3(X(:,1), X(:,2), X(:,3), 25, U(:,3), 'filled');
view(-20,20)
axis equal; axis off; axis vis3d;
q = 100;
for it=1:q
    view((it-1)/q*360,20); drawnow;
    % mysaveas(it,'eig2');
end


if n<=600
    [I,J] = find(K);
    clf; hold on;
    for k=1:length(I)
        i = I(k); j = J(k);
        plot3( X([i j],1), X([i j],2), X([i j],3), 'k' );
    end
    scatter3(X(:,1), X(:,2), X(:,3), 25, t(:,1), 'filled');    
    axis equal;
    axis off;
    view(-20,20);
    for it=1:q
        view((it-1)/q*360,20); drawnow;
        mysaveas(it,'wire');
    end
end

