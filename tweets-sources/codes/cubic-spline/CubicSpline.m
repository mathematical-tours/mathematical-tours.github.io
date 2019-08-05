addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


% sampling points and values
n = 12;
x = rand(n,1);
y = rand(n,1);

x = cumsum(.5 + rescale(rand(n,1), 0, 1));
x = rescale(x,.05,.95);
y = rescale(y, .1,.9);

% grid
p = 1024;
g = (0:p-1)/p;
I = floor(x*p);

% min |f(I)-y|^2 + lambda*<Delta*f,f>
S = sparse(I,I,ones(n,1),p,p);
Y = zeros(p,1); Y(I)=y;

% Laplacian in 1D
B = [-ones(p,1),2*ones(p,1),-ones(p,1)];
Delta = spdiags(B,[-1 0 1],p,p);
Delta(1,end) = -1;
Delta(end,1) = -1;


% smoothing spline
r = 1; lmax = 50;
r = 2; lmax = 100000;
r = 3; lmax = 100*500000;
q = 50;
lambda_list = linspace(.01,lmax,q);

for it=1:q    
    lambda = lambda_list(it);
    f = (S+lambda*Delta^r)\Y;
    %
    s = (it-1)/(q-1);
    clf; hold on;
    plot(g, f, 'LineWidth', 2, 'Color', [s 0 1-s]);
    stem(I/p,y, 'k', 'filled', 'MarkerSize', 6);
    axis([0 1 0 1])
    box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end

% AutoCrop(rep, 'anim-'); 