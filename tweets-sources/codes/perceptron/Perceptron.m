%%
% Test for approximation of 1D function using a perceptron with a single
% hidden layer.

addpath('../toolbox/');
rep = MkResRep();

% sigmoid
phi = @(x)1./(1+exp(-x));
phiD = @(x)phi(x).*(1-phi(x));
% f(x,a,b,c) = sum_i c_i phi(a_i*x-b_i) 
% min_{a,b,c} E(a,b,c) = 

n = 256; % #samples
p = 10; % #hidden layer
x = linspace(-1,1,n)';

% target function 
y = sin(pi*x);

Phi  = @(a,b)phi(x*a'-ones(n,1)*b');
PhiD = @(a,b)phiD(x*a'-ones(n,1)*b');

% loss
E = @(a,b,c)1/2*norm(Phi(a,b)*c-y)^2;
% grad
R = @(a,b,c)Phi(a,b)*c-y;
GEc = @(a,b,c)Phi(a,b)' * R(a,b,c);
GEa = @(a,b,c)( PhiD(a,b)' .* ( c * R(a,b,c)' ) ) * x;
GEb = @(a,b,c)( PhiD(a,b)' .* ( c * R(a,b,c)' ) ) * ones(n,1);

% init
a0 = 10+randn(p,1);
c0 = randn(p,1);
b0 = 2*rand(p,1)-1;

a = a0; b = b0; c = c0;
tau = .1;
for i=1:niter
    % display
    clf; hold on;
    plot(x, Phi(a,b)*c, 'LineWidth', 2);
    plot(x, y, 'k', 'LineWidth', 2);
    axis([-1 1 -1.5,1.5]); box on; set(gca, 'XTick', [], 'YTick', []);
    % descent
end
