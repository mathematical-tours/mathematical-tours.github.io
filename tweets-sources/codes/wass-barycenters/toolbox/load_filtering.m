function [blur,K] = load_filtering(filt_method, n)

% load_filtering - load a Gaussian filtering operator callback
%
%   blur = load_filtering(filt_method, n);
%
%   Copyright (c) 2014 Gabriel Peyre

% blurring
cconv = @(a,Fb)real(ifft2( fft2(a).*Fb ));
% symmetric extension
ext = @(x)[x; x(end:-1:1,:)];
ext = @(x)ext(ext(x)')';
rest = @(x)x(1:end/2,1:end/2);

% 1D Gaussian kernel
if nargout>1
    [Y,X] = meshgrid(1:n,1:n);
    K = @(mu)exp( -(X-Y).^2/(2*mu^2) );
end

% metric options -- for fft filtering only
gaussian_opt.exponent = 2; % put =2 for W2, =1 for W1
gaussian_opt.threshold = inf;

switch filt_method
    case 'fftper'   %% here assuming 2D 
        blur = @(x,mu)cconv( x, fft2(gaussian(n,mu,gaussian_opt)) );
    case 'fftsym'   %% here assuming 2D 
        blur = @(x,mu)rest( cconv( ext(x), fft2(gaussian(2*n,mu,gaussian_opt)) ) );
    case 'imgaussian'
        blur = @(x,mu)imgaussian(x,mu,mu*50);
    case 'gaussianiir'
        numsteps = 8;
        blur = @(x,mu)gaussianiir2d(x,mu, numsteps);
    case 'matrix'         %% only for 1D
        blur = @(x,mu)K(mu)*x;
    otherwise
        error('Unknown filtering method');
end

% check for symmetry of the operator
if 0
dotp = @(x,y)sum(x(:).*y(:));
a = randn(n,n); b = randn(n,n);
mu = 10;
v = dotp(blur(a,mu),b); 
w = dotp(a,blur(b,mu));
err = abs(v-w)/abs(v);
if err>1e-5
    warning('The operator is non symetric');
end
end