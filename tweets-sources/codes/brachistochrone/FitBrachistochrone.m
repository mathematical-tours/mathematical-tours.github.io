% display of a Brachistochrone

addpath('../toolbox/');
rep = MkResRep();

t = linspace(-pi,pi,1024*2);
Br = @(R)R*(t+sin(t)+pi) - 1i*R*(1+cos(t));

xmax = 2;
q = 10;
Rlist = linspace(.32,.5,q);
clf; hold on;
X = [];
for it=1:q
    r = (it-1)/(q-1);
    B = Br(Rlist(it));
    [~,i] = min(abs(real(B)-xmax));
    plot(B(1:i), 'Color', [r 0 1-r], 'LineWidth', 2);
    plot(B(i), '.', 'Color', [r 0 1-r], 'MarkerSize', 25);
    X(it) = imag(B(i));
end
plot(0,0, '.k', 'MarkerSize', 25);
axis equal; axis([0 xmax*1.01 -1 0]); axis off;
saveas(gcf, [rep 'curves-brachi.eps'], 'epsc');


x = linspace(0,xmax,100);
clf; hold on;
for it=1:q
    r = (it-1)/(q-1);
    plot(x, x/xmax*X(it), 'Color', [r 0 1-r], 'LineWidth', 2);
    plot(xmax, X(it), '.', 'Color', [r 0 1-r], 'MarkerSize', 25);
end
plot(0,0, '.k', 'MarkerSize', 25);
axis equal; axis([0 xmax*1.01 -1 0]); axis off;
saveas(gcf, [rep 'curves-lines.eps'], 'epsc');