function [f1,E,z] = perform_tv_denoising(f,lambda, options)

% perform_tv_denoising - solve the ROF TV denoising problem
%
% compute the solution f1 of 
%   min_f1 1/2*|f-f1|^2 + lambda*|grad(f1)|_1
% It uses a projected gradient descent to solve
%   min 1/2*|div(z) + f/lambda|^2  s.t. |z|_inf<=1
% and output
%   f1 = f+lambda*div(z)
%
% Copyright (c) 2016 Gabriel Peyre

mynorm = @(x)norm(x(:));

N = size(f,1);
[grad,div,lapl] = load_grad(N);
L = 16; % norm of laplacian

options.null = 0;
niter = getoptions(options, 'niter', 1000);

ampl = @(z)sqrt(sum(z.^2,3));
Proj1 = @(a)min(a,1) ./ max(a,1e-10);
Proj = @(z)z .* repmat( Proj1(ampl(z)), [1 1 4] );
tau = 1.8/L;
tau = 1.8/L;

z = getoptions(options, 'z0', []);
if isempty(z)
    z = zeros(N,N,4);
end

E = [];
for i=1:niter
    progressbar(i,niter);
    u = div(z) + f/lambda;
    E(i) = mynorm(u);
    z = Proj( z + tau*grad( u ) );
end
% denoising result.
f1 = f+lambda*div(z);

end
