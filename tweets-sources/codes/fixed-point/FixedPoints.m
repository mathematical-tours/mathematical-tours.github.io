%%
% Display 1D fixed point iterations.

addpath('../toolbox/');
rep = MkResRep();

eta = - .4;
eta = -.31;
eta = -.5;
eta = .32;

f = @(x)3*x .* (1-x);
f = @(x)x + eta *sin(x*1*pi*2);

x = linspace(-.3,1.3,1024);

m = 20;

clf;
subplot(2,1,1); hold on;
plot(x, x);
plot(x, f(x));
subplot(2,1,2);
plot(x, abs( (f(x)-f(x+1e-5))/(1e-5) ));

q = 50;
ulist = linspace(-.1,1.1,q);
for it = 1:q
    u = ulist(it);    
    clf; hold on;
    plot(x,x, 'k', 'LineWidth', 1);
    plot(x,f(x), 'k', 'LineWidth', 2);
    for i=1:m
        c = (i-1)/(m-1);
        v = f(u);
        plot([u u],[u v], 'color',[c 0 1-c], 'LineWidth', 1);
        plot([u v],[v v], 'color',[c 0 1-c], 'LineWidth', 1);
        u = v;
    end
    box on;
    axis square;  axis([min(x) max(x) min(x) max(x)]); 
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end