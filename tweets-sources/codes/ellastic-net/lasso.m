function [W,lambda_list,lmax] = lasso(X0,y0, lambda_list,rho, options)

% R(x) = lambda*( (1-rho)*|x|_1 + rho/2*|x|^2 )
%      = s*|x|_1 + r/2*|x|^2 
%   s=(1-rho)*lambda, r=rho*lambda

lmax = max(X0'*y0);

options.null = 0;
lambda_scale = getoptions(options, 'lambda_scale', .1);


% q = 400;
% lambda_list = getoptions(options, 'lambda_list', lmax*linspace(lambda_scale,1e-3,q) );
q = length(lambda_list);

C = X0'*X0;
u = X0'*y0;

Soft = @(x,s)max(abs(x)-s,0).*sign(x);
% with ellastic net
%    min 1/2*|x-y|^2 + s*|x| + r/2*|x|^2
% (r+1)/2*x^2 - xy + ... =
%   min 1/2*|x-y/(r+1)|^2 + s/(r+1)*|x|
Ellas = @(y,s,r)Soft(y/(r+1),s/(r+1));

ISTA = @(w,s,r,tau)Ellas( w-tau*( C*w-u ), s*tau, r*tau );

tau = 1.5/norm(X0)^2;
p = size(C,1);
% regul
W = []; E = [];
w = zeros(p,1);
niter = getoptions(options, 'niter', 1000);
for iq=1:q
    progressbar(iq,q);
    lambda = lambda_list(iq);
    %
    s = (1-rho)*lambda; r=rho*lambda;
    % ISTA %
    for i=1:niter
        w = ISTA(w,s,r,tau);
    end
    W(:,iq) = w; % bookkeeping
    % E(iq) = Delta(X0*w,y1);
end


end