%%
% Convergence of shot noise to Gaussian process.

name = 'iso';
name = 'aniso';

addpath('../toolbox/');
rep = MkResRep(name);


% image size
n = 512; 

% shot noise pattern
x = [0:n/2,-n/2+1:-1]'/n;
[Y,X] = meshgrid(x,x);

switch name
    case 'iso'
        s = .03;
        h = exp(-(X.^2+Y.^2)/(2*s^2));
        vmax = 15;
    case 'aniso'        
        s = .02;
        f = 20;
        h = exp(-(3*X.^2+ 1.5*Y.^2)/(2*s^2)) .* cos(2*pi*X*f);
        vmax = 12;
        
end

q = 75;
plist = round( linspace(1,20000) );

convol = @(x,y)real(ifft2(fft2(x).*fft2(y)));

for i=1:q
    p = plist(i);
    A = zeros(n);
    rand('state', 123); randn('state', 123);
    I = ceil(rand(p,1)*n*n);
    A(I) = randn(p,1);
    % convolution
    f = convol(A,h);
    clf;
    imagesc(f);
    caxis([-1 1]*vmax);
    colormap gray(256);
    axis image; axis off; 
    drawnow;
end