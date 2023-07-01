%  test for nystrom kernel approximation, in 1D.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


p = 2; 
s = .03;
k = @(x,y)exp(-abs( (x(:)-y(:)')/s ).^p);

n = 1000; % for display
m = 100; % samples
x = linspace(0,1,n)';
x0 = rand(m,1);
% x0 = linspace(0,1,m)';

for s=2:m
    z = x0(1:s);
    K = k(x,z) * inv(k(z,z)) * k(z,x);
    clf;
    imagesc(K);
    drawnow;
end


return;

% data to interpolate
n = 12;
rand('state', 12);
randn('state', 12);
x = rand(n,1);

x = cumsum(1+rand(n,1));
x = rescale(x,.05,.95);
y = rand(n,1);

lambda = .001;


r = 0;
v = k(z,x) * pinv( k(x,x) + lambda*eye(n) ) * y;
clf; hold on;
plot(z,v, 'color', [r 0 1-r], 'LineWidth', 2);
plot(x,y, 'k.', 'MarkerSize', 25);
axis([0 1 -.2 1.2]); box on;
set(gca, 'PlotBoxAspectRatio', [1 2/3 1],'XTick', [], 'YTick', []);
drawnow;