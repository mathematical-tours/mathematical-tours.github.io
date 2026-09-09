[im, cm] = imread('maurice.png');
colormap(cm);

[N,P] = size(im);

im = double(im);


NN = 64;
x_start = 80;
y_start = 80;
im_extr = im(x_start:x_start+NN-1, y_start:y_start+NN-1);

% image d'origine
SUBPLOT(2,2,1);
image(im);
axis image;
axis off;
title('(a) Image d''origine');
% un petit carré
h = line( [x_start x_start], [y_start y_start+NN-1] );
set( h, 'LineStyle', ':' );
h = line( [x_start+NN-1 x_start+NN-1], [y_start y_start+NN-1] );
set( h, 'LineStyle', ':' );
h = line( [x_start x_start+NN-1], [y_start y_start] );
set( h, 'LineStyle', ':' );
h = line( [x_start x_start+NN-1], [y_start+NN-1 y_start+NN-1] );
set( h, 'LineStyle', ':' );

% extrait
SUBPLOT(2,2,2);
image(im_extr);
axis image;
axis off;
title('(b) Image extraite');

% calcule la correlation
correl = xcorr2(im, im_extr);
m = max(max(correl));
correl = correl*256/m;
SUBPLOT(2,2,3);
image(1:N,1:N,correl(1:N,1:N));
axis image;
axis off;
title('(c) Correlation');

% renormalisation
% kernel = ones(NN,NN);
% nor = xcorr2(im.^2, kernel);
% correl = correl./sqrt(nor);

correl = normxcorr2(im_extr, im);

m = max(max(correl));
correl = correl*256/m;
SUBPLOT(2,2,4);
image(1:N,1:N,correl(1:N,1:N));
axis image;
axis off;
title('(d) Correlation normalisée');

correl = correl(1:N,1:N);

% recherche du max
[Q,QQ] = size(correl);
m = 0;
xs = 0;
ys = 0;
for i=1:N
for j=1:N
	if correl(i,j)>m
		xs = i;
		ys = j;
		m = correl(i,j);
    end
end
end

h = line( [xs xs], [0 400] );
set( h, 'LineStyle', ':' );
h = line( [0 400], [ys ys] );
set( h, 'LineStyle', ':' );

saveas(gcf, '../correlation-2d', 'eps')
saveas(gcf, '../correlation-2d', 'png')