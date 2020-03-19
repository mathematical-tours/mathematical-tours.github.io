% iterative convolution
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

h = [1 1];
h = [1 29];
h = h/sum(h);

q = 50;
x = 1;
for n=1:q
    s = (n-1)/(q-1);
    a = (n-1)/2;
    clf;
    col = [s 0 1-s];
    u = (-a:a)/sqrt(n); m = sum(u.*x);
    bar( u-m, x, 'EdgeColor', col, 'FaceColor', col);
    set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
    axis([-2 2 0 .2]); % max(x)*1.02
    drawnow;
    mysaveas(n);
    %
    x = conv(x, h);
    
end