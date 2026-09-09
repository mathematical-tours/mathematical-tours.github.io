clear;
N = 512;
x = 2*pi*(0:N-1)/(N-1);
z = exp(i*x);
r1 = 1./(z-1);
r2 = 1/2*(z+1)./(z-1);
r3 = 1/2*(1+4*z+z.^2)./(z.^2-1);

subplot(1,3,1);
plot(x, abs(r1), 'k');
axis([0,2*pi, 0, 2]);
axis square;
title('H^{(1)}(\omega)');

subplot(1,3,2);
plot(x, abs(r2), 'k');
axis([0,2*pi, 0, 2]);
axis square;
title('H^{(2)}(\omega)');

subplot(1,3,3);
plot(x, abs(r3), 'k');
axis([0,2*pi, 0, 2]);
axis square;
title('H^{(3)}(\omega)');

saveas(gcf, '../meth-quadratures-filtre-recursif', 'eps')
saveas(gcf, '../meth-quadratures-filtre-recursif', 'png')