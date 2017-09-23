%%
% Test for 2D FFT.


name = 'hibiscus';
F = imread([name '.png']);

F = sum(F,3)/255;
clf; imagesc(F); axis image; colormap gray;

G = fft2(F);
clf;
imagesc(log(abs(G)+1e-10));

clf;
imagesc(fftshift( log(abs(G)+1e-10) ));

%% 
% Same with windowing

n = size(F,1);
t = (0:n-1)'/n;
h = sin(pi*t).^2;

clf; plot(h); axis tight;

Fw = (h*h').*F;
clf; imagesc(Fw);


Gw = fft2(Fw);
clf;
imagesc(fftshift( log(abs(Gw)+1e-10) ));
