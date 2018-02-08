%%
% Test for Metropolis-Hastings for a discrete distribution.

% dimension
n = 100;
% #samples
q = 1000000;

niter = 500;

% initial locations
X = ceil(rand(q,1)*n);

% target distribution
mu = rand(n,1);
mu = [1:n/2, n/2:-1:1]';
mu = mu/sum(mu);

% radius for move
r = 2;

for i=1:niter
    nu = hist(X,1:n); nu = nu/sum(nu);
    clf; hold on;
%    bar(1:n, mu, 1, 'r', 'EdgeColor', 'r');
%    bar(1:n, nu, 'b');
    area(1:n, nu);
    plot(1:n, mu, 'r', 'LineWidth', 2);
    axis tight; drawnow;
    % proposal move
    V = ceil(rand(q,1)*(2*r+1))-r-1;
    X1 = mod(X + V-1,n)+1;
    % accept/reject
    R = mu(X1) ./ mu(X);
    U = rand(q,1)>R;
    X = X1.*(1-U) + X.*U;
end
