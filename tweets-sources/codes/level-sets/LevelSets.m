%%
% Display of levelsets.

rep = '../results/level-sets/';
[~,~] = mkdir(rep);

addpath('../toolbox/');

name = 'hibiscus';
name = 'lisa';
name = 'pacman1';
n = 512;

% gaussian blur
t = [0:n/2,-n/2+1:-1]';
G = @(s)exp(-t.^2/(2*s^2)); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));


f0 = rescale(sum(load_image(name, n),3));
if f0(1,1)>f0(end/2,end/2)   
    f0 = 1-f0;
end
f = GFilt(f0,17);
f = rescale(f);
% equalize
% [fS,I] = sort(f(:)); f(I) = linspace(0,1,length(f(:)));


Quant = @(x,q)min(floor( rescale(x)*q  ), q-1);

% display quantized colormap
t = linspace(0,1,n);
for r=[5 10 20]
    clf; hold on;
    imagesc(t,t,f);
    M = linspace(0,1,r); M([1 end]) = [];
    contour(t,t,f,M, 'r');
    colormap(parula(r-1));
    caxis([0 1]); axis ij;
    axis image; axis off;
    saveas(gcf, [rep name '-' num2str(r) '.png'], 'png');
    % 3D plot
    clf; 
    surf(Quant(f,r-1));
    shading interp;
    view(-30,60); axis off;
    camlight;
    saveas(gcf, [rep name '-' num2str(r) '-3d.png'], 'png');
end