%%
% Test for betti number from simplicial homology

addpath('../toolbox/');
rep = MkResRep('betti');
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


n = 200;  % #points
n0 = 512;
name = 'shape2';
f = rescale(sum(load_image(name, n0),3));
if f(1)==1
    f = 1-f;
end
I = find(f(:)==1);
I = I(randperm(length(I)));
[x,y] = ind2sub([n0 n0],I(1:n));
x = (x-1)/(n0-1);
y = (y-1)/(n0-1);

X = [x,y]';
D = distmat(X,X);

s = .05;
q = 40;
slist = linspace(0,.13,q);

Beta = []; 

for it=1:q
    progressbar(it,q);
    t = (it-1)/(q-1);
    col = [t 0 1-t];
    s = slist(it);
    %
    [i,j] = find(sparse(D<s));
    E = [i(i<j), j(i<j)];
    F = []; 
    F = TriangleFromEdges(E);
    T = []; 
    % T = TriangleFromEdges(F);
    % compute betti number
    [B,beta] = ComputeBoundaries(E,F,T);
    Beta(:,end+1) = beta;
    %
    % do the display
    clf; hold on;
    patch('Faces',F,'Vertices',X','FaceColor', .1*col+.9, 'EdgeColor', col);
    plot([x(E(:,1)) x(E(:,2))]', [y(E(:,1)) y(E(:,2))]', 'color', col, 'LineWidth', 2);
    plot(x,y, 'k.', 'MarkerSize', 20);
    axis equal; axis([0 1 0 1]); axis off;
    drawnow;
    mysaveas('evol',it);
end

Beta(:,1) = Beta(:,2);
for it=1:q
clf; hold on;
plot(slist, Beta(1,:), 'b', 'LineWidth', 2);
plot(slist(it), min(Beta(1,it),40), 'b.', 'MarkerSize', 40);
axis([0 max(slist), 0 40]); box on;
set(gca, 'FontSize', 15, 'PlotBoxAspectRatio', [1 1/2 1]);
drawnow;
mysaveas('beta0',it);

clf; hold on;
plot(slist, Beta(2,:), 'r', 'LineWidth', 2);
plot(slist(it), Beta(2,it), 'r.', 'MarkerSize', 40);
axis tight; box on;
set(gca, 'FontSize', 15, 'PlotBoxAspectRatio', [1 1/2 1]);
mysaveas('beta1',it);
end