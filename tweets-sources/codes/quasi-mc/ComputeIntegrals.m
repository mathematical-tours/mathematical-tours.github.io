addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

d = 3;
p = sobolset(3);

n = 20000;
p = scramble(p,'MatousekAffineOwen');
x = net(p,n);
x = 2*x-1;
y = 2*rand(n,d)-1;

v0 = 4/3*pi;

vx = 8*cumsum(sum(x.^2,2)<1) ./ (1:n)';
vy = 8*cumsum(sum(y.^2,2)<1) ./ (1:n)';

q = 150;
for it=1:q
    m = round(it/q*n);
    clf; 
    semilogy((abs(vx-v0)), 'LineWidth', 1, 'color', .5*[0 0 1] + .5);
    hold on;
    semilogy((abs(vy-v0)), 'LineWidth', 1, 'color', .5*[1 0 0] + .5);
    semilogy((abs(vx(1:m)-v0)), 'b', 'LineWidth', 2);
    semilogy((abs(vy(1:m)-v0)), 'r', 'LineWidth', 2);
    semilogy(m, (abs(vx(m)-v0)), 'b.', 'MarkerSize', 25);
    semilogy(m, (abs(vy(m)-v0)), 'r.', 'MarkerSize', 25);
    axis([0 n 1e-3 1]);
    % axis([0 n -3 0]);
    set(gca, 'XTick', [], 'YTick', [1e-3 1e-2 1e-1], 'FontSize', 13);
    box on;
    drawnow;
    mysaveas(it);
end


return;

it = q;
for it=1:q
    
end