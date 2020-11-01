% shapiro inequality


addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

t = linspace(.000001,3,200);
[y,x] = meshgrid(t,t);
A = 1./(x+y) + x./(1+y) + y./(1+x);


A = 1./(x+y) + x./(1+y) + y./(1+x);

u = 1; z = 1; 
A = x./(y+z+u) + y./(z+u+x) + z./(u+x+y) + u./(x+y+z);

A = rescale(A);
% display quantized colormap
r = 15; % #levellines
clf; hold on;
imagesc(t,t,A');
contour(t,t,A',linspace(0,1,r), 'k');
colormap(parula(r-1));
caxis([0 1]);
axis image; axis off;