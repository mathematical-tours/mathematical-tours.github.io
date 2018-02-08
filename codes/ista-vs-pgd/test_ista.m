lambda = 1;
tau = 1;
N = 201;
x = linspace(-2.5,2.5,N)';
P  = @(x)max(x,0);
Pm = @(x)max(-x,0);
eta = randn(N,1);

Soft = @(x,k)sign(x).*P(abs(x)-k);

R = P(P(x)-tau*(eta+lambda)) - ...
    P(Pm(x)-tau*(-eta+lambda));

clf; hold on;
plot(R);
plot(Soft(x-tau*eta,lambda*tau), 'r.:');


clf; hold on;
plot(R - Soft(x-tau*eta,lambda*tau), 'r');


%%
% implement lasso

N = 200; P = round(N/2);
A = randn(P,N);
y = randn(P,1);
lambda = .1;
tau = 1/norm(A)^2;
x = zeros(N,1);
niter = 300;
E_fista = [];
for i=1:niter
    eta = A'*(A*x-y);
    x = Soft(x-tau*eta,lambda*tau);
    if mod(i,10)==0
    clf;
    plot(x); drawnow;
    end
    E_fista(i) = 1/2*norm(y-A*x)^2+lambda*norm(x,1);
end


%%
% implement PGD

Proj  = @(x)max(x,0);
u = zeros(N,1);
v = zeros(N,1);
Coefs = [];
E_pgd = [];
% tau = 1/(2*norm(A)^2);
for i=1:niter
    eta = A'*(A*(u-v)-y);
    %
    u = Proj(u-tau*(eta+lambda));
    v = Proj(v-tau*(-eta+lambda));
    Coefs(i) = sum(u~=0 & v~=0);
    %    
    x1 = u-v;
    if mod(i,10)==0
    clf; hold on;
    plot(u-v)
    plot(u, 'g.');
    plot(-v, 'r.'); 
    drawnow;
    end
    E_pgd(i) = 1/2*norm(y-A*(u-v))^2+lambda*norm(u-v,1);
end

clf; hold on;
plot([E_fista(:),E_pgd(:)]);
axis tight;
legend('Fist','PGD');

clf;
plot(( E_fista-E_pgd ) / max(E_fista));
axis tight;