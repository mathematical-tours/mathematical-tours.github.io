%%
% Test for Ribs complex

addpath('../toolbox/');
rep = MkResRep('ribs');
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);

X = [x,y]';
D = distmat(X,X);

s = .05;
q = 100;
slist = linspace(0,.3,q);

for it=1:q
    t = (it-1)/(q-1);
    col = [t 0 1-t];
    %
    s = slist(it);
    [i,j] = find(sparse(D<s));
    E = [i(i<j), j(i<j)];
    F = TriangleFromEdges(E);
    %
    clf; hold on;
    patch('Faces',F,'Vertices',X','FaceColor', .1*col+.9, 'EdgeColor', col);
    plot([x(E(:,1)) x(E(:,2))]', [y(E(:,1)) y(E(:,2))]', 'color', col, 'LineWidth', 2);
    plot(x,y, 'k.', 'MarkerSize', 20);
    axis equal; axis([0 1 0 1]); axis off;
    drawnow;
    mysaveas('evol',it); 
end




