[im, cm] = imread('maurice_noise.png');
colormap(cm);

im_svg = im;

% image d'origine
SUBPLOT(1,2,1);
image(im);
axis image;
axis off;
title('Image bruitée');

% on étend la taille des matrices pour faire une
% convolution
[n,p] = size(im);
nn = 3;	% taille du support du filtre
t = n+nn; % nouvelle taille des matrices
% on agrandi pour calculer la convolution
im(t, t) = 0;
f = zeros(t, t);
f(1,1) = 4;
f(2,1) = 2;
f(1,2) = 2;
f(t,1) = 2;
f(1,t) = 2;
f(3,1) = 1;
f(2,2) = 1;
f(1,3) = 1;
f(t-1,1) = 1;
f(t,2) = 1;
f(1,t-1) = 1;
f(2,t) = 1;
f(t,t) = 1;
f(4,1) = 0.5;
f(3,2) = 0.5;
f(2,3) = 0.5;
f(1,4) = 0.5;
f(t-2,1) = 0.5;
f(t-1,2) = 0.5;
f(t,3) = 0.5;
f(1,t-2) = 0.5;
f(2,t-1) = 0.5;
f(3,t) = 0.5;
f(t,t-1) = 0.5;
f(t-1,t) = 0.5;
m = sum(sum(f));
f = f./m;

im = real( ifft2(fft2(im).*fft2(f)) );
im = im(2:n,1:n);
SUBPLOT(1,2,2);
image(im);
axis image;
axis off;
title('Après filtrage par une gaussienne');

saveas(gcf, '../../images/filtre-2d', 'eps');
saveas(gcf, '../../images/filtre-2d', 'png');
