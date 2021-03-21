addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);

n = 25;
rand('state', 12);
randn('state', 12);
x = rand(n,1);
y = sin(x*2*pi)+randn(n,1)*.15;


p = 1;
U = @(x)(x-1).*(x>.5) + (x+1).*(x<-.5) + x.*(x>=-.5).*(x<=.5);
d = @(x,z)abs(U(x(:)-z(:)'));
K = @(x,z,s)exp( -d(x,z).^p/(s^p) );


clf;
plot(x,y, 'k.', 'MarkerSize', 25);

z = linspace(0,1,1000);

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
    mysaveas('anim',it);
    
    clf;
    plot(z, K(1/2,z,s), 'color', [r 0 1-r], 'LineWidth', 2);
    axis([0 1 -.05 1.05]); box on;
    set(gca, 'PlotBoxAspectRatio', [1 1/3 1], 'XTick', [], 'YTick', []);
    mysaveas('kernel',it);
end

