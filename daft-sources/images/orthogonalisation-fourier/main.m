% on veut que (f,f(.+1),f(.+2)...) soit une famille orthogonale

global g;

n = 20;
% f = rand(n,1);
c = [(n/2+1):n, 1:n/2]';
x = (-n/2+1):n/2;

f = zeros(n,1);
f(1) = 1;
f(2) = 0.5;
g = orthog(f);
subplot(2,2,1);
plot(x, f(c), 'k*:');
title('Fonction originale');
axis tight;
subplot(2,2,2);
plot(x, g(c), 'k*:');
title('Fonction orthogonalisée');
axis tight;


f = zeros(n,1);
f(1) = 1;
f(2) = 0.9;
g = orthog(f);
subplot(2,2,3);
plot(x, f(c), 'k*:');
title('Fonction originale');
axis tight;
subplot(2,2,4);
plot(x, g(c), 'k*:');
title('Fonction orthogonalisée');
axis tight;


saveas(gcf, '../orthogonalisation-fourier', 'eps')
saveas(gcf, '../orthogonalisation-fourier', 'png')

f = exp(2*i*pi/n* (0:n-1) )';