%%
% Fourier cristalographics

addpath('../toolbox/');
rep = MkResRep();

n = 256;
M = zeros(n);
t = linspace(-1,1,n);

[Y,X] = meshgrid(t,t);

K = 7;
Z = exp(2i*pi/K*(1:K)');

sizeI = [1024,1024];
spacing = 20;
nPts = 0;
showIter = 0;
U = poissonDisc(sizeI,spacing,nPts,showIter);
Z = ( U(:,1) + 1i*U(:,2) )/sizeI(1);


% click selection
Z = [];
while true
    clf; hold on;
    imagesc(t,t,real(M)'); 
    colormap parula(256);
    plot(real(Z), imag(Z), 'k.', 'MarkerSize', 25);
    axis([-1 1 -1 1]);
    axis image; axis off;
    [a,b,button] = ginput(1);
    if button==3
        break;
    end
    Z(end+1) = a+1i*b;
    % compute Fourier transform
    rho = 5;
    M = zeros(n);
    for s=1:length(Z)
        M = M + exp( rho * 2i*pi*( X*real(Z(s)) + Y*imag(Z(s)) ) );
    end
    imageplot(abs(M));
    
end