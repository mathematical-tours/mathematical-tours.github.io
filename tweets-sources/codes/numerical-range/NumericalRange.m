%%
% Test for the  numerical range


addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep  'anim-' znum2str(it,3) '.png']);

n = 12;

% random matrix
A = randn(n);
% normal matrix
[U,~] = qr(randn(n));
B = U * diag(randn(n,1)+1i*randn(n,1)) * U';

q = 50;
for it=1:q
    t = (it-1)/(q-1);
    At = (1-t)*A+t*B;
    clf; hold on;
    [RERANGE, IMRANGE, NRADIUS] = nrange(At,400, [t 0 1-t]);
    plot(eig(At), '.', 'color', [t 0 1-t], 'MarkerSize', 30);
    axis equal; axis off;
    drawnow;
    mysaveas(it);
end