%%
% Compute exactly the Voronoi cells.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

n = 3;
x = rand(n,1)+1i*rand(n,1);

col = rand(3,500);

q=150;

for it=1:q
    A = [0; 1; 1+1i; 1i]; % add dummy vertices
    xx = [x; (A -.5-.5i)*500];
    [V,C] = voronoin([real(xx), imag(xx)]);
    B = zeros(length(x),1);
    clf; hold on;
    for i=1:length(x)
        I = C{i}([1:end,1]);
        % intersection with bounding box
        p1 = polyshape([0 0 1 1],[1 0 0 1]);
        p2 = polyshape(V(I,1),V(I,2));
        pp = intersect(p1,p2);
        W = pp.Vertices;
        % compute barycenter
        [u,v] = centroid(pp);
        B(i) = u+1i*v;
        %
        % plot( W([1:end 1],1), W([1:end 1],2), 'color', rand(3,1), 'LineWidth', 2);
        fill( W([1:end 1],1), W([1:end 1],2), col(:,i)' );
    end
    plot(x, 'r.',  'MarkerSize', 10);
    axis equal; axis([0 1 0 1]);
    axis on; box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(it);
    x = B;
    % add points
    if 1 % mod(it,2)==1
        K = 1; 
        x(end+1:end+K) = .5+.5i + (randn(K,1)+1i*randn(K,1))*.05;
    end
end

return;
V = VX + 1i*VY;



% Display Delaunay
clf; hold on;