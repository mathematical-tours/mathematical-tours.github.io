
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


n = 40;

M = randn(n); M = (M+M')/sqrt(2);
M = .5*M/sqrt(n);
% hist(eig(M), 100);

u = randn(n,1);

tmax = .05;

q = 200;
niter = 4000;
ndisp = round(1 + (niter-1)*linspace(0,1,q).^5); 
ndisp = unique(ndisp);
clist = linspace(0,1,niter);

k = 1;
clf; hold on;
for i=1:niter
    t = (i-1)/(niter-1);
    z = eig( M + tmax * 1i*t*u*u' );
    plot(real(z), imag(z), '.', 'MarkerSize', 10, 'color', [clist(i) 0 1-clist(i)]);
    set(gca, 'PlotBoxAspectRatio', [1 1 1], 'XTick', [], 'YTick', []);
    axis([-1,1,0,.3]); box on;
    if i==ndisp(k)
        drawnow;
        k = k+1;
        mysaveas(k);
    end
end


