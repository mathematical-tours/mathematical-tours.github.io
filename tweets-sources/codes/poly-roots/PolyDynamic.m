%%
% Display roots of polynomials


rep = '../results/poly-roots/';
[~,~] = mkdir(rep);
if not(exist('test'))
    test=1;
end

% p5, ... p 0
p = [1 0 0 0 0 -1];

q = 201;
u = 1.6;
x = linspace(-u,u,q);
[Y,X] = meshgrid(x,x);
Z = X+1i*Y;

F = polyval(p,Z);
R = roots(p);

clf; hold on;
imagesc(x,x,abs(F));
contour(x,x,sqrt(abs(F)), 25, 'k', 'LineWidth', 2);
plot(imag(R),real(R), 'r.', 'MarkerSize', 20);
colormap parula(256);
axis off;

% radius
r = [0.5 0.5 0.5 0.5 0.5 1];
% center
c = [1 0.1 0.1i -.2 0.1+.1i 0];


% radius
r = [0 0 2 0 0 0];
c = [1 0 0 0 0 1];

m = 100;
nls = 9;
for i=1:m
    t = (i-1)/m;
    p = c + exp(2i*pi*t)*r;
    %   
    F = polyval(p,Z);
    R = roots(p);
    %
    clf; hold on;
    imagesc(x,x,abs(F));
    contour(x,x,(abs(F)), linspace(0,4,nls), 'k', 'LineWidth', 2);
    plot(imag(R),real(R), 'r.', 'MarkerSize', 20);
    colormap(parula(nls-1)); caxis([0 4]);
    axis([-u u -u u]);
    axis off;
    drawnow;
    saveas(gcf, [rep 'anim-' num2str(i) '.png'], 'png');
end
