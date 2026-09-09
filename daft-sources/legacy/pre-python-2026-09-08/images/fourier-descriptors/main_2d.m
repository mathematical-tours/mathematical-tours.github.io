% test of trigonometric interpolation of polygons
clf;
[x,y] = ginput;
n = length(x);
z = x+1i*y;

p = 1025;
fz = fft(z);
q = ceil(n/2);
zz = [fz(1:q); zeros(p-n,1); fz((q+1):n)];
zz = p/n*ifft(zz);
z1 = zz;


fz = fft(x);
zz = [fz(1:q); zeros(p-n,1); fz((q+1):n)];
xx = p/n*real( ifft(zz) );

fz = fft(y);
zz = [fz(1:q); zeros(p-n,1); fz((q+1):n)];
yy = p/n*real( ifft(zz) );
zz = xx+1i*yy;
z2 = zz;

z = [z; z(1)];
zz = [zz; zz(1)];
plot(real(z), imag(z), 'k*-', real(zz), imag(zz), 'k:');
axis tight;

