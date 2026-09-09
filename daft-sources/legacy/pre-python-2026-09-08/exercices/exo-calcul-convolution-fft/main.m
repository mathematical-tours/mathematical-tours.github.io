% test de la convolution optimisée.

n = 128;
m = 8;
f = rand(n,1);
g = rand(m,1);

y = convol(f,g);
yy = conv(f,g); % la fonction matlab
disp( norm(y-yy) );