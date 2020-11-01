% test de le FHT. Comparaison avec la FFT.

n = 16;
p = 512;    % taille avec zero pad.
f = [(0:n/2-1)';(n/2:-1:1)'];

f = [f;zeros(p-n,1)];
g = fht(f);
g = [g(p/2+1:p);g(1:p/2)];
h = real(fft(f));
h = [h(p/2+1:p);h(1:p/2)];


subplot(1,2,1);
plot(0:n-1, f(1:n), 'k*:');
title('Signal original');
axis square;
axis tight;

subplot(1,2,2);
x = (-p/2+1:p/2)*(n/p);
plot(x, g, 'k-', x, h, 'k-.');
legend('FHT', 'FFT');
axis square;
axis tight;

title('Transformées de Fourier et de Hartley');
saveas(gcf, '../../images/comparaison-fht-fft', 'eps');
saveas(gcf, '../../images/comparaison-fht-fft', 'png');