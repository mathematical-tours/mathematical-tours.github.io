n = 15;
nn = (n+1)/2;
p = 512;

f = rand(n,1);

ff = fft(f);
ff = [ff(1:nn); zeros(p-n,1); ff(n-nn+2:n)];
ff = real(ifft(ff))*(p/n);

hh = fht(f);
hh = [hh(1:nn); zeros(p-n,1); hh(n-nn+2:n)];
hh = fht(hh)*(1/n);


h = 1/p;
% xx = (0:h:1-h)*(n-1);
xx = (0:(p-1))*(n/p);
x = 0:n-1;

plot(x,f,'k*', xx,ff,'k', xx,hh,'k-.');
axis([0,n-1,-0.1,1.1]);
legend('Signal original', 'Interpolation Fourier', 'Interpolation Hartley');

% on obtient la meme interpolation
norm(ff-hh)

% saveas(gcf, '../interpolation-hartley', 'eps');
% saveas(gcf, '../interpolation-hartley', 'png');

