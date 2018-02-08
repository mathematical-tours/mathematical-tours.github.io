%%
% Test for Metropolis-Hastings for a discrete distribution.

% dimension
n = 10;
% #samples
q = 1000;

niter = 100;

% initial locations
X = ceil(rand(q,1)*n);

% target distribution
mu = rand(n,1); mu = mu/sum(mu);

% radius for move
r = 1;

for i=1:niter
    clf; hold on;
    bar(1:n, mu, 'r'); axis tight;
    hist(X,1:n); axis tight;
    % proposal move
    V = ceil(rand(q,1)*(2*r+1))-r-1;
    X1 = mod(X + V,n)+1;
    % accept/reject
    R = mu(X1) ./ mu(X);
    U = rand(q,1)<R;
    X = X1.*(1-U) + X.*U;
end
