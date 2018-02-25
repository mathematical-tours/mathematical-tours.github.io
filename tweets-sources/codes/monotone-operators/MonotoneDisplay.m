%%
% Display of monotone operators

rep = '../results/monotone-operators/';
[~,~] = mkdir(rep);
addpath('../toolbox/');

r = 15; % #levellines

% display background
q = 501;
T = linspace(-1,1,q);
[Y,X] = meshgrid(T,T);
% display vectors
m = 11;
t = linspace(-.85,.85,m);
[y,x] = meshgrid(t,t);

% sym
F = X.^2+Y.^2; % function
U = x; V = y;
clf; hold on;
imagesc(T,T,F');
contour(T,T,F',linspace(min(F(:)), max(F(:)),r), 'k');
quiver(x,y,U,V, 'r', 'LineWidth', 2);
colormap(parula(r-1));
caxis([min(F(:)) max(F(:))]);
axis image; axis off;
saveas(gcf, [rep 'sym.png'], 'png');

% skew
F = X.*Y; % function
U = y; V = -x;
clf; hold on;
imagesc(T,T,F');
contour(T,T,F',linspace(min(F(:)), max(F(:)),r), 'k');
quiver(x,y,U,V, 'r', 'LineWidth', 2);
colormap(parula(r-1));
caxis([min(F(:)) max(F(:))]);
axis image; axis off;
saveas(gcf, [rep 'skew.png'], 'png');

