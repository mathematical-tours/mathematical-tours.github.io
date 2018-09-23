function PlotMST(X)


D = distmat(X,X);
% compute delaunay triangulation graph.
DT = delaunayTriangulation(X');
T = DT.ConnectivityList;
% build graph
[MST, pred] = graphminspantree(sparse(D));
[I,J] = find(MST);
% display graph
clf; hold on;
triplot( T, X(1,:), X(2,:), '-', 'Color', [1 1 1]*0, 'LineWidth', 1, 'MarkerSize', 30 );
triplot( T, X(1,:), X(2,:), 'r.', 'LineWidth', 1, 'MarkerSize', 25 );
for k=1:length(I)
    plot(X(1,[I(k) J(k)]), X(2,[I(k) J(k)]), 'r-', 'LineWidth', 3);
end
axis equal; axis([0 1 0 1]);
% axis off;
axis on; box on; set(gca,'Xtick', [],'Ytick', []);


end