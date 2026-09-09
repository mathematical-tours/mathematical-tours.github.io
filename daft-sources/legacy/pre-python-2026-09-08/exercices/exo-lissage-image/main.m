[im, cm] = imread('maurice.png');
n = length(im); s = 0.01;
f = calcul_filtre(n,s);
m = max(max(f));
y = filter2(f,im);
image(y); colormap(cm);
axis off; axis image;