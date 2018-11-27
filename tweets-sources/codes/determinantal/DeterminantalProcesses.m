%%
% test for translation invariant sampling of determinantal processes


addpath('../toolbox/');
rep = MkResRep();

n = 256*2;
t = [0:n/2,-n/2+1:-1]'/n;
[Y,X] = meshgrid(t,t);


sx = .2; sy = .03;
h = exp(-((X/sx).^2+(Y/sy).^2)/2);

sx = .2; sy = .01;
h = exp(-((X/sx).^2+(Y/sy).^2)/2) + exp(-((X/sy).^2+(Y/sx).^2)/2);

s = .06*3;
h = exp(-(X.^2+Y.^2)/(2*s^2));


FC = abs(fft2(h));
FC = FC/max(FC(:));

rand('state', 123);
randn('state', 123);
X = DPixP2DNew(n,n,FC);