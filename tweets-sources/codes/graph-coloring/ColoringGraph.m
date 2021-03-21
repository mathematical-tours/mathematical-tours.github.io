
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

n = 150;
Z = rand(n,1)+1i*rand(n,1);
q = 70;

for it=1:q
    n = it+2;
    z = Z(1:n);
    % delaunay graph
    DT = delaunayTriangulation([real(z), imag(z)]);
    I = [DT.ConnectivityList(:,1);DT.ConnectivityList(:,2);DT.ConnectivityList(:,3)];
    J = [DT.ConnectivityList(:,2);DT.ConnectivityList(:,3);DT.ConnectivityList(:,1)];
    A = sparse(I(:),J(:),ones(length(I(:)),1));
    A = max(A,A');
    [I,J,~] = find(A);
    %
    G = graph(I,J);
    C = GraphColoringJohnson(G)+1;
    
    col = distinguishable_colors(25); col(4,:)=[];
    clf; hold on;
    for s=1:length(I)
        i = I(s); j = J(s);
        plot(z([i j]), 'k-', 'LineWidth', 2);
    end
    for k=1:max(C(:))
        S = find(C==k);
        plot(z(S), '.', 'color', col(k,:), 'MarkerSize', 30);
    end
    axis equal; axis([0 1 0 1]); axis off;
    drawnow;
    mysaveas(it);
end