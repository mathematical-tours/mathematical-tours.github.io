%%
% Display variational formulation of |.|=min_eta


addpath('../toolbox/');
rep = MkResRep();

q = 50;
elist = linspace(.8,.01,q);
elist = .001 + .8 * linspace(1,0,q).^4;


x = linspace(-1,1,1024);

for it=1:q
    clf; hold on;
    for jt=1:it
        eta = elist(jt);
        s = (jt-1)/(q-1);
        plot( x, 1/2*( x.^2/eta + eta ), 'Color', [s 0 1-s]*.3+.7, 'LineWidth', 2  );
    end
    plot(x,abs(x), 'k', 'LineWidth', 2);
    plot( x, 1/2*( x.^2/eta + eta ), 'Color', [s 0 1-s], 'LineWidth', 3  );
    plot( [-1 1]*eta, [1 1]*eta, '.', 'MarkerSize', 25, 'Color', [s 0 1-s], 'LineWidth', 2  );
    box on;
    set(gca, 'XTick', [], 'YTick', []);
    axis equal; axis([-1 1 -.02 1.5]);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end

x = linspace(-1,1,1024);
eta = linspace(1e-3,1,512);
[Eta,X] = meshgrid(eta,x);
F = 1/2*(X.^2./Eta+Eta);

r = 15; % #levellines
vmax = 1.5;
clf; hold on;
imagesc(x,eta,F');
contour(x,eta,F',linspace(0,vmax,r), 'k');
plot(x, abs(x), 'r', 'LineWidth', 2);
colormap(parula(r-1));
caxis([0 vmax]);
axis image; axis off;
saveas(gcf, [rep 'eta-approx.png']);

