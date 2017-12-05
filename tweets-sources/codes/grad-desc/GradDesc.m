%%
% Test for gradient descent of a simple 2D function.

n = 200;
xt = linspace(-.3,1.5,n);
yt = linspace(-.3,1.5,n);
[Y,X] = meshgrid(yt,xt);

F = Rosenb(X,Y);

v = sort(F(:));
u = v(round(1:n*n/20:end));

clf; hold on;
imagesc(xt,yt,F');
contour(xt,yt,F', u, 'k');
colormap parula(256);
axis equal; box on;



% initial point
x = 1.2; y = 0.8;
niter = 100;
U = []; V = [];
F1 = [];
tau = 1e-4;
for i=1:niter
    U(end+1) = x; V(end+1) = y;
    [f,gx,gy] = Rosenb(x,y);
    F1(end+1) = f;
    x = x - tau*gx; 
    y = y - tau*gy;
end


clf; hold on;
imagesc(xt,yt,F');
contour(xt,yt,F', 15, 'k');
plot(U,V, 'r.-', 'LineWidth', 2, 'MarkerSize', 23);
colormap parula(256);
axis image; axis off;
drawnow;


% saveas(gcf, [rep 'gd-' num2str(a) '.png'], 'png');