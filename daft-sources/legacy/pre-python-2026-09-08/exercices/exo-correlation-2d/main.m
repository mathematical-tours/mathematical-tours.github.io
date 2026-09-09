[im, cm] = imread('maurice.png');
colormap(cm);
im = double(im);
N = length(im);

P = 64;
x_start = 80;
y_start = 80;
im_extr = im(x_start:x_start+P-1, y_start:y_start+P-1);

% image d'origine
SUBPLOT(2,2,1);
image(im);
axis image;
axis off;
title('(a) Image d''origine');

% extrait
SUBPLOT(2,2,2);
image(im_extr);
axis image;
axis off;
title('(b) Image extraite');



correl = normxcorr2(im_extr, im);
m = max(max(correl));
correl = correl*256/m;
SUBPLOT(2,2,3);
image(1:N,1:N,correl(1:N,1:N));
axis image;
axis off;
title('(d) Correlation normalisée');

correl = correlation_normalisee(im,im_extr);
m = max(max(correl));
correl = correl*256/m;
SUBPLOT(2,2,4);
image(1:N,1:N,correl);
axis image;
axis off;
title('(d'') Correlation normalisée');