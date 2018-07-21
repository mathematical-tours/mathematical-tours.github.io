addpath('../toolbox/');
rep = MkResRep('dynamic');

n = 1000;
U = (2*rand(n)-1)*sqrt(2);
U = randn(n);

nmin = 1;
nmax = 80;
q = 80;
nlist = round(linspace(nmin,nmax,q));

for i=1:q
    t = (i-1)/(q-1);
    k = nlist(i);
    clf; hold on;
    plot(eig(U(1:k,1:k)/sqrt(k)), '.', 'MarkerSize', 25, 'color', [t 0 1-t]);
    axis equal; 
    axis([-1 1 -1 1]*1.05); 
    axis off;
    plot( exp(2i*pi*linspace(0,1,100)), 'k' );
    drawnow;
    saveas(gcf, [rep 'wigner-' znum2str(i,2), '.png']);
end

% AutoCrop(rep, ['wigner-']);