% Roger W. Brockett
% Dynamical systems that sort lists, diagonalize matrices, and solve linear programming problems. Linear Alg. Appl., 146:79?01,
% 1991

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep  'anim-' znum2str(it,3) '.png']);

com = @(X,Y)X*Y-Y*X;


n = 80;
X0 = randn(n); X0 = X0+X0';

[Q,R] = qr(randn(n));
X0 = Q*diag(linspace(-1,1,n))*Q';

q = 130;
dt = .1*3;
dt = .1+2*linspace(0,1,q).^2;
A = diag(linspace(0,1,n));
X = X0;
for it=1:q
    X = X + dt(it) * com(X,com(X,A));
    clf;
    % imagesc(X); 
    imagesc(log(.1+abs(X)));
    axis image; axis off;
    drawnow;
    mysaveas(it);
end