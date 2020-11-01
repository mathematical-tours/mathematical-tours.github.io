
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

t = linspace(-1,1,200);
[Y,X] = meshgrid(t,t);

p = .5;
normp = @(x,y)( abs(x).^p + abs(y).^p ).^(1/p);
if p<1
normp = @(x,y)( abs(x).^p + abs(y).^p );
end
q = 75;

for it=1:q
    s = (it-1)/q;
    x0 = .8*cos(2*pi*s); y0 = .8*sin(2*pi*s);
    A = normp(X+x0,Y+y0) ./ ( normp(X,Y) + normp(x0,y0) );
    clf;
    r = 15; % #levellines
    clf; hold on;
    imagesc(t,t,A');
    cmax = .99;
    contour(t,t,A',linspace(0,cmax,r), 'k');
    colormap(parula(r-1));
    caxis([0 cmax]);
    plot(x0,y0,'r.', 'MarkerSize', 25);
    axis image; axis off;
    drawnow;
    mysaveas(it);
end