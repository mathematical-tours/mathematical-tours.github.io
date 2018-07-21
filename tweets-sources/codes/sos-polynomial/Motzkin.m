%% Motzkin polynomial

addpath('../toolbox/');
rep = MkResRep();

fmax = 1;

a = 1.35;
t = linspace(-a,a,511);
[Y,X] = meshgrid(t,t);
f = X.^4 .* Y.^2 + X.^2 .* Y.^4 - 3*Y.^2.*X.^2 + 1;

% display quantized colormap
r = 15; % #levellines
clf; hold on;
imagesc(t,t,f');
contour(t,t,f',linspace(0,fmax,r), 'k');
colormap(parula(r-1));
caxis([0 fmax]);
axis image; axis off;
saveas(gcf, [rep 'motzkin-2d.png'], 'png');


% display
clf; hold on;
surf(t,t,f');
shading interp;
% view(-65,40);
view(3)
colormap(parula(r-1));
caxis([0,fmax]);
axis([-a a -a a 0 fmax]);
camlight; axis off;
saveas(gcf, [rep 'motzkin-3d.png'], 'png');