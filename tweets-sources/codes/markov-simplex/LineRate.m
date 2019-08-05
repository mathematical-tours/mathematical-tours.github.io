%%
% Test for linear rate for Perron Frobenius fixed points.

rep = '../results/hilbert-metric/';
[~,~] = mkdir(rep);

% dimension
n = 40; 

u = ones(n,1);

name = 'rand';
name = 'test1';
name = 'test2';
name = 'isotropic';

normalize = @(x)x/sum(x);

% input matrix
switch name 
    case 'rand'
        K = rand(n)+20*eye(n);
    case 'isotropic'
        K = ones(n)+20*eye(n);
end

% make it stochastic, K'*u=1
c = randn(n);

epsilon = 1;
K = exp(-c/epsilon);


% K = rand(n) + eye(n)*1000;

K = K*diag(1./(K'*u));


% hilber metric
%   log max_{ij} (xi/xj)/(yi/yj)

% H = @(x,y) log(max(max( (x*(1./x')) ./ (y*(1./y')) ) ) );
Var = @(x)max(x(:))-min(x(:)); % variation norm
H = @(x,y)Var(log(x)-log(y));

% invariant proba
[V,S] = eig(K); 
[~,I] = sort(diag(S)); 
v = V(:,I(end)); v=v/sum(v);


x = normalize(rand(n)); % each column is a point of the simplex
E = [];
niter = 100;
for i=1:niter
    E(i) = H(x,v);
    x = K*x;
end

[lambda,eta] = ContractionRatio(K);

[lambda,eta] = ContractionRatio(K); lambda
[lambda,eta] = ContractionRatio(diag(rand(n,1))*K*diag(rand(n,1))); lambda

% other way to compute
S = -Inf; 
for i=1:n
    for j=i:n
        S = max(S,H(K(:,i),K(:,j)));
    end
end
eta1 = exp(S);

tau = @(eta)( sqrt(eta)-1 ) / ( sqrt(eta)+1 );
lambda1 = tau(eta1)

clf; hold on;
plot(log10(E/E(1)));
plot(1:niter, log10(lambda)*(1:niter), '--' );
axis tight; 

