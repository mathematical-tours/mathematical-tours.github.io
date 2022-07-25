addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);

dotp = @(a,b)sum(a(:).*b(:));

tol = 1e-10;
niter = 100;


n = 200;
n = 30;


A = diag( (1:n) );
b = A*ones(n,1);



A = diag( (1:n).^2 );
A = diag( (1:n) );
b = ones(n,1);


p = 3;
A = diag( ( .1 + (0:n-1) ).^p );
b = A*ones(n,1);

niter = 60;
x = zeros(n,1);
r = b-A*x;
p = r;
E = [];
for it=1:niter
    r1 = r;
    E(it) = norm(r(:));
    Ap = A*p;
    alpha = dotp(r,r) / dotp(p,Ap);
    x = x + alpha*p;
    r = r-alpha*Ap;
    beta = dotp(r,r)/dotp(r1,r1);
    p = r + beta*p;  
    
    x1 = x;
        warning off;
    if 1
       % [x1,~,~] = cgs(A,b,1e-20,it);  
      % [x1,~,~] = gmres(A,b,[],1e-20,it);  
       [x1,~,~] = minres(A,b,1e-20,it);  
    end
     warning on;


    s = (it-1)/(niter-1);
    clf;
    bar(1:n,x1, 'FaceColor', [s 0 1-s]);
    axis([.5 n+.5 0 1.15]);
    set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
    drawnow;
%    mysaveas('anim',it)
end
