function [u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(mu,nu,c,epsilon,options)

% sinkhorn_log - stabilized sinkhorn over log domain with acceleration
%
%   [u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(mu,nu,c,epsilon,options);
%
%   mu and nu are marginals.
%   c is cost
%   epsilon is regularization
%   coupling is 
%       gamma = exp( (-c+u*ones(1,N(2))+ones(N(1),1)*v')/epsilon );
%
%   options.niter is the number of iterations.
%   options.tau is an avering step. 
%       - tau=0 is usual sinkhorn
%       - tau<0 produces extrapolation and can usually accelerate.
%
%   options.rho controls the amount of mass variation. Large value of rho
%   impose strong constraint on mass conservation. rho=Inf (default)
%   corresponds to the usual OT (balanced setting). 
%
%   Copyright (c) 2016 Gabriel Peyre

options.null = 0;
niter = getoptions(options, 'niter', 1000);
tau  = getoptions(options, 'tau', -.5);
verb = getoptions(options, 'verb', 1);
rho = getoptions(options, 'rho', Inf);
tol = getoptions(options, 'tol', 0);

lambda = rho/(rho+epsilon);
if rho==Inf
    lambda=1;
end

N = [size(mu,1) size(nu,1)];
H1 = ones(N(1),1);
H2 = ones(N(2),1);

ave = @(tau, u,u1)tau*u+(1-tau)*u1;


lse = @(A)log(sum(exp(A),2));
M = @(u,v)(-c+u*H2'+H1*v')/epsilon;

% kullback divergence
H = @(p)-sum( p(:).*(log(p(:)+1e-20)-1) );
KL  = @(h,p)sum( h(:).*log( h(:)./p(:) ) - h(:)+p(:) );
KLd = @(u,p)sum( p(:).*(exp(-u(:))-1) );
dotp = @(x,y)sum(x(:).*y(:));

err = [];
u = zeros(N(1),1); 
v = zeros(N(2),1);
gamma = exp(M(u,v));
Wprimal = []; Wdual = [];
for i=1:niter
    if verb==1
        progressbar(i,niter);
    end
    u1 = u;
    u = ave(tau, u, ...
        lambda*epsilon*log(mu) - lambda*epsilon*lse( M(u,v) ) + lambda*u );
    v = ave(tau, v, ...
        lambda*epsilon*log(nu) - lambda*epsilon*lse( M(u,v)' ) + lambda*v );
    % coupling 
    gamma = exp(M(u,v));
    if rho==Inf % marginal violation
        Wprimal(i) = dotp(c,gamma) - epsilon*H(gamma);
        Wdual(i) = dotp(u,mu) + dotp(v,nu) ...
            - epsilon*sum( gamma(:) );
        err(i,1) = norm( sum(gamma,2)-mu );
    else % difference with previous iterate
        Wprimal(i) = dotp(c,gamma) - epsilon*H(gamma) ...
            + rho*KL(sum(gamma,2),mu) ...
            + rho*KL(sum(gamma,1),nu);
        Wdual(i) = -rho*KLd(u/rho,mu) - rho*KLd(v/rho,nu) ...
            - epsilon*sum( gamma(:) );
        err(i,1) = norm(u(:)-u1(:), 1);
    end
    if err(i,1)<=tol
        break
    end
end

end