%%
% Evolution of a voronoi cell.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

name = 'shapes';

n0 = 512;
n = 150;
f0 = rescale(-load_image([name], n0))>.5; f0 = f0(:,:,1);
f = rescale(-load_image([name], n)) >.5; f = f(:,:,1);

t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
I = find(f(:));

g = X(:) + Y(:)*1i; 
z = X(I) + Y(I)*1i; % point of the object


x = .5 + 1i*.5;

% distant to object
dz = reshape( min(abs(g - transpose(z)), [], 2), [n n]);
dx = reshape( abs(g - x), [n n]);
V = (dx<=dz);

clf; hold on;
