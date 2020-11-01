[im, cm] = imread('maurice.png');
n = length(im);
p = ceil(n/2*(sqrt(2)-1));
x = ones(n+2*p,n+2*p)*255;
sel = (p+1):(n+p);
x(sel,sel) = im;
y = fft_rotation(x,0.2);
