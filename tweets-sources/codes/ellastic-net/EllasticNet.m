%%
% Display Ellastic Net Path

rep = '../results/ellastic-net/';
[~,~] = mkdir(rep);

addpath('../toolbox/');

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


% generate data
n = 10; 
p = round(n);
randn('state', 123);
A = randn(p,n);
y = randn(p,1);
% normalize columns
A = A ./ repmat( sqrt(sum(A.^2,1)), [p 1] );

lmax = max(abs(A'*y))*1.2;
q = 400;
lambda_list = linspace(.001,1,q)*lmax;

rho = 0; % lasso
rho = 1; % ridge
options.niter = 4000;
[W] = lasso(A,y, lambda_list,rho); % (A,y, options);

clf; hold on;
plot(lambda_list, W', 'LineWidth', 2);
plot([0 lmax], [0 0], 'k--', 'LineWidth', 2);
a = max(abs(W(:)));
axis([0 lmax -a a]); box on;
SetAR(2/3);
set(gca, 'XTick', [], 'YTick', []);
saveas(gcf, [rep 'path' num2str(round(rho*10)) '.eps'], 'epsc');
