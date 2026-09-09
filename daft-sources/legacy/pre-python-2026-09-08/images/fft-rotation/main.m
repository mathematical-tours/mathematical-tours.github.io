[im, cm] = imread('maurice.png');
n = length(im);
p = ceil(n/2*(sqrt(2)-1));
x = ones(n+2*p,n+2*p)*255;
x((p+1):(n+p),(p+1):(n+p)) = im;
nbr = 6; rot = pi/(4*(nbr-1)); 
for r = 0:nbr-1
    y = fft_rotation(x,r*rot);
    subplot(1,nbr,r+1);
    image(y); colormap(cm);
    axis off; axis image;
end
saveas(gcf, '../../images/fft-rotation', 'eps')
saveas(gcf, '../../images/fft-rotation', 'png')
