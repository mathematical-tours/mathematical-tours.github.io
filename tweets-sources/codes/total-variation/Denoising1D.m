n = 1024; %% increasing this slow down the convergence
 
x = linspace(0,1,n)';
f = abs(x-.5)<=.3;
 
grad = @(x)x-x([2:end 1]);
div  = @(x)x([end 1:end-1])-x;
L = 4; % norm of laplacian
 
%projected gradient descent on the dual
Proj = @(z)z ./ max(abs(z),1);
tau = 1.8/L;
 
lambda = 10;
 
mfun = @(x) lambda*norm(grad(x),1)+ norm(x-f)^2/2;
 
z = zeros(n,1);
niter = 1000;
E = [];
fval = [];
for i=1:niter
    u = div(z) + f/lambda;
    E(i) = norm(u, 'fro')^2;
    z = Proj( z + tau*grad( u ) );
    % denoising result.
    f1 = f+lambda*div(z);
    if mod(i,50)==1
        clf; hold on;
        plot(x, f, 'k--');
        plot(x, f1, 'r', 'LineWidth', 2);
        axis([0 1 0 1]);
        drawnow;
    end
    fval(i) = mfun(f1);
end
 
 
 
%% VarPro
 
%sparse gradient matrix
e = ones(n,1);
D =  spdiags([ e -1*e ],0:1,n,n);
D(n,1)=-1;
 
test = randn(n,1);
disp(['should be zero:', num2str(norm(D*test-grad(test)))])
 
 
v = 0.01*ones(n,1);
fval_v = [];
for i=1:niter
    
    a = (lambda*(D*D')+spdiags(v.^2,0,n,n))\(D*f);
    g = v - a.^2.*v;
    v = v - tau*g;
    
    f1v = -lambda*D'*a + f;
    
    if 0 % mod(i,2)==1
        clf; hold on;
        plot(x, f, 'k--');
        plot(x, f1v, 'r', 'LineWidth', 2);
        axis([0 1 0 1]);
        drawnow;
    end
   
    fval_v(i) = mfun(f1v);
end
a = (lambda*(D*D')+spdiags(v.^2,0,n,n))\(D*f);
f1v = -lambda*D'*a + f;
 
%% plot objective
figure(3)
clf
objmin = min(min(fval_v), min(fval));
loglog(fval - objmin)
hold on
loglog(fval_v - objmin, 'r')
legend('PCG', 'VarPro')