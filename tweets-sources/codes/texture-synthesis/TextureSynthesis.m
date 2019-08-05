% 'filename': the image file containing the sample image (the texture to grow)
% 'winsize': the edge length of the window to match at each iteration (the window is (winsize x winsize) )
% (newRows, newCols): the size of the output image


addpath('../toolbox/');
rep = MkResRep();

m = 64;
n = 256;


name = 'cheetah';
name = 'plants';

f = imread([name '.jpg']);
f = f(1:m,1:m,:);

imwrite(rescale(double(f)), [rep 'original.png' ]);

w = 7;
[g, h] = synth(f, w, n, n, rep);
