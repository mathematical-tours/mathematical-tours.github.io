% test for group lassp trajectories

addpath('../toolbox/');
rep = MkResRep();

name = 'group';
name = 'ridge';
name = 'lasso';

n = 10; p = 5;
randn('state', 1234);
A = randn(n,p)+1i*randn(n,p);
x = zeros(p,1);
if strcmp(name, 'lasso')
    A = [real(A)-imag(A), 1i*(real(A)+imag(A))];
    x = zeros(2*p,1);
end

y = randn(n,1)+1i*randn(n,1);
lambdam = max(abs(A'*y));

nb = 100;
lambda_list = linspace(lambdam,0,nb);

niter = 5000;
Xa = [];
for i=1:nb
    lambda = lambda_list(i);
    progressbar(i,nb);
    switch name
        case 'group'
            tau = 1.2/norm(A)^2;
            X = ista(A,y,x,lambda,tau, niter);
            x1 = X(:,end);
        case 'lasso'
            tau = 1.2/norm(A)^2;
            X = ista(A,y,x,lambda,tau, niter);
            x = X(:,end);
            x1 = x(1:end/2) + 1i*x(end/2+1:end);
        case 'ridge'
            x1 = (A'*A + lambda*eye(p)) \ (A'*y);
    end
    Xa(:,end+1) = x1;
end

clf;
plot3(lambda_list, real(Xa'), imag(Xa'), 'LineWidth', 2);