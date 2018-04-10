%%
% Display of function on torus.


rep = '../results/mean-median/';
[~,~] = mkdir(rep);

n = 256;

t = linspace(0,1,256);
[V,U] = meshgrid(t,t);
V = V;
U = U;


R = 1; r = .5;

Torus = @(U,V)deal((R+r*cos(2*pi*V)).*cos(2*pi*U), ...
    (R+r*cos(2*pi*V)).*sin(2*pi*U), ...
    sin(2*pi*V));

[Y,Y,Z] = Torus(U,V);

r = 16; % #levellines
% display 3D torus

name = 'xy';
switch name
    case 'y'
        f = Y;
    case 'xy'
        f = X.*Y;
end

[cm,c]=contour(t,t,f,linspace(min(f(:)),max(f(:)),r));

clf; hold on;
surf(Z,X,Y, f);
view(-80,10); axis equal;
shading interp;
camlight; axis off;
%
while not(isempty(cm))
    k = cm(2,1); cm(:,1) = [];
    u = cm(1,1:k); v = cm(2,1:k); cm(:,1:k) = [];
    [x,y,z] = Torus(v,u);
    plot3(z,x,y,'k', 'LineWidth', 3);
end
colormap(parula(r-1));
caxis([min(f(:)),max(f(:))]);
saveas(gcf, [rep 'torus-' name '.png'], 'png');

%% 
% 2-D

r = 24;
s=2.5;
f = peaks(s*(2*U-1),s*(2*V-1));

% display quantized colormap
clf; hold on;
imagesc(t,t,f);
contour(t,t,f,linspace(min(f(:)),max(f(:)),r), 'k', 'LineWidth', 2);
colormap(parula(r-1));
caxis([min(f(:)),max(f(:))]);
axis image; axis off;
saveas(gcf, [rep 'peaks-2d.png'], 'png');

