addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

p = sobolset(2);

n = 600;
p = scramble(p,'MatousekAffineOwen')
x = net(p,n);
% x = rand(n,2);

q = 150;
for it=1:q
    m = round(n*it/q);
    clf;
    scatter(x(1:m,1), x(1:m,2), 15, (1:m)/n, 'filled'); 
    box on;
    set(gca, 'XTick', [], 'YTick', []);
    axis equal;
    axis([0 1 0 1]);
    drawnow;
    % mysaveas(it);
end