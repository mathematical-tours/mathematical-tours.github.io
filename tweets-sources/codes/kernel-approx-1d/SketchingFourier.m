%%
% Use random Fourier features.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);

n = 12;
rand('state', 12);
randn('state', 12);
x = rand(n,1);

x = cumsum(1+rand(n,1));
x = rescale(x,.05,.95);
y = rand(n,1);
% y = sin(x*2*pi)+randn(n,1)*.15;


U = @(x)(x-1).*(x>.5) + (x+1).*(x<-.5) + x.*(x>=-.5).*(x<=.5);
U = @(x)x;
d = @(x,z)abs(U(x(:)-z(:)'));

s = .05;
K = @(x,z)exp( -d(x,z).^2/(2*s^2) );
lambda = .01;

clf;
plot(x,y, 'k.', 'MarkerSize', 25);

z = linspace(0,1,1000);


r = 0;
v = K(z,x) * pinv( K(x,x) + lambda*eye(n) ) * y;
clf; hold on;
plot(z,v, 'color', [r 0 1-r], 'LineWidth', 2);
plot(x,y, 'k.', 'MarkerSize', 25);
axis([0 1 -.2 1.2]); box on;
set(gca, 'PlotBoxAspectRatio', [1 2/3 1],'XTick', [], 'YTick', []);
drawnow;

q = 2000; % #Fourier frequencies

m = 70;
qmax = 2000;
qlist = round( 1+(qmax-1)*linspace(0,1,m).^2 );
qlist = unique(qlist);
m = length(qlist);

omega = randn(qmax,1)*1/s;
omega(1) = 0; % be sure to sample the constant
phi = @(z,q)exp(1i*omega(1:q)*z(:)');
Ka = @(x,z,q)phi(x,q)'*phi(z,q)/q;

for it=1:m    
    q = qlist(it);
    r = (it-1)/(m-1);
    va = real( Ka(z,x,q) * pinv( Ka(x,x,q) + lambda*eye(n) ) * y );
    clf; hold on;
    plot(z,v, 'k:', 'LineWidth', 2);
    plot(z,va, 'color', [r 0 1-r], 'LineWidth', 2);
    plot(x,y, 'k.', 'MarkerSize', 25);
    axis([0 1 -.2 1.2]); box on;
    set(gca, 'PlotBoxAspectRatio', [1 2/3 1],'XTick', [], 'YTick', []);
    drawnow;
    mysaveas('approx',it);
    %
    g = K(z,0.5);
    ga = real( Ka(z,0.5,q) );
    % g = g/max(g);
    clf; hold on;
    plot(z, g, 'k:', 'LineWidth', 2);
    plot(z, ga, 'color', [r 0 1-r], 'LineWidth', 2);
    axis([0 1 -.1 1.1]); box on;
    set(gca, 'PlotBoxAspectRatio', [1 1/3 1], 'XTick', [], 'YTick', []);
    mysaveas('kernel',it);

end



return;

q = 70;
slist = linspace(.01,.2,q);
lambda = .01;
for it=1:q
    r = (it-1)/(q-1);
    s = slist(it);
    v = K(z,x,s) * pinv( K(x,x,s) + lambda*eye(n) ) * y;
    clf; hold on;
    plot(z,v, 'color', [r 0 1-r], 'LineWidth', 2);
    plot(x,y, 'k.', 'MarkerSize', 25);
    axis([0 1 -1.3 1.3]); box on;
    set(gca, 'PlotBoxAspectRatio', [1 2/3 1],'XTick', [], 'YTick', []);
    drawnow;
    % mysaveas('anim',it);
    
    clf;
    plot(z, K(1/2,z,s), 'color', [r 0 1-r], 'LineWidth', 2);
    axis([0 1 -.05 1.05]); box on;
    set(gca, 'PlotBoxAspectRatio', [1 1/3 1], 'XTick', [], 'YTick', []);
    %mysaveas('kernel',it);
end