rep = 'results/';
[~,~] = mkdir(rep);

a = 1.5;
n = 256;
t = linspace(-a,a,n);
[Y,X] = meshgrid(t,t);

name = 'ellpsoid';
name = 'paraboloid';
name = 'quadric';
name = 'quadric-ncvx';

switch name
    case 'ellpsoid'
        P = X.^2 + 2*Y.^2; % -.5*Y).^2+1.5*Y.^2;
    case 'paraboloid'
        P = 2*X.^2-Y.^2;
    case 'quadric'
        P = X.^4+Y.^4;
    case 'quadric-ncvx'
        P = X.^4+Y.^4-4*X.^2.*Y.^2;
end

lw = 2;
q = 20;
clf; hold on;
imagesc(t,t,P);
contour(t,t,P,q,'k', 'LineWidth', lw);
axis equal; axis([-a a -a a]);
axis off;
colormap parula(256);
saveas(gcf, [rep name '.png']);