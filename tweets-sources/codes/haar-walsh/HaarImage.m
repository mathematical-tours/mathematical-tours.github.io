%%
% Haat of an image.

addpath('../toolbox/');
rep = MkResRep('image');

name = 'hibiscus';

n = 256;
f = load_image(name,n);
f = rescale(sum(f,3));
imwrite(clamp(f), [rep name '-original.png']);

myscale = @(x).5 + .15*x/std(x(:));

f = [(f(1:2:end,:)+f(2:2:end,:))/2;(f(1:2:end,:)-f(2:2:end,:))/2];


g = f; 
u = 1:n/2; v = n/2+1:n;
g(v,:) = myscale(g(v,:));
imwrite(clamp(g), [rep name '-w1.png']);

f = [(f(:,1:2:end)+f(:,2:2:end))/2,(f(:,1:2:end)-f(:,2:2:end))/2];
g = f; 

g(u,v) = myscale(g(u,v));
g(v,v) = myscale(g(v,v));
g(v,u) = myscale(g(v,u));

imwrite(clamp(g), [rep name '-w2.png']);