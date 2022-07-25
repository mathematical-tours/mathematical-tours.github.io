addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

n = 2048*4+1;
x = linspace(-.7,.7,n);

delta = (max(x)-min(x))/n;
dx  = @(x)1/delta*(x([2:end 1])-x);
dxS = @(x)1/delta*(x-x([end 1:end-1]));
epsilon = 1e-25;
abs1 = @(x)sqrt(x.^2 + epsilon^2);
lapl = @(x,p)dxS( dx(x) .* abs1(dx(x)).^(p-2) );

s = .2;
f = exp(-x.^2/(2*s^2));

clf;
plot(x,f, 'k', 'LineWidth',2);
axis([min(x) max(x) -.05 1.05]);
box on; set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
drawnow;
saveas(gcf, [rep 'function.png']);

q = 100;
plist = linspace(.5,3,q);

g = lapl(f,max(plist));
vmin = min(g)*1.05;
vmax = max(g)*1.05;

for it=1:q
    clf; hold on;
    for jt=1:it
        s = (jt-1)/(q-1);
        col = [s 0 1-s];
        p = plist(jt);
        g = lapl(f,p);
        plot(x(2:end-1),g(2:end-1), 'color', .5 + .5*col, 'LineWidth', 1);
    end    
    plot(x(2:end-1),g(2:end-1), 'color', col, 'LineWidth', 3);
    axis([min(x) max(x) vmin vmax]);
    box on;
    set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(it);
end

