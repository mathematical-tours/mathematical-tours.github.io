%% 
% code for density estimator.

% create save repertory
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


n = 400;
gauss = @(n,m,s)m+s*(randn(n,1)+1i*randn(n,1));
z = [gauss(n/2,.6+.6i,.07); gauss(n/2,.2+.3i,.1)];
z = rand(n,1)+1i*rand(n,1);


clf;
plot(z, 'r.', 'MarkerSize', 10);
axis equal;
axis([0 1 0 1]);

p = 1024;
q = 70;
slist = linspace(0.001,.2,q);
slist = 0.01 + .15*linspace(0,1,q).^2;
for it=1:q
    s = slist(it)*p;
    h = parzen2d(real(z),imag(z),p,s);
    h = max(h,0);
    h = h/max(h(:));
    % display
    r = 15; % #levellines
    t = linspace(0,1,p);
    clf; hold on;
    imagesc(t,t,h');
    contour(t,t,h',linspace(0,1,r), 'k');
    colormap(parula(r-1));
    caxis([0 1]);
    plot(z, '.', 'MarkerSize', 5, 'color', [1 .5 .5]);
    axis image; axis off;
    drawnow;
    mysaveas(it);
end