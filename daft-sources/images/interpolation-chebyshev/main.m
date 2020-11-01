N = 11;
alpha = 0.3;
NN = 1000;

% fonction à interpoler
xx = -1:2/(NN-1):1;
yy = 1./(alpha^2+xx.^2);

% coeff du polynôme de la lagrange
x = -1:2/(N-1):1;
y = 1./(alpha^2+x.^2);
[P,S] = polyfit(x,y,N-1);
yy_lagrange = polyval(P,xx);

% coeff du polynôme de la tchebyshev
x = 0:(N-1);
x = -cos( x*pi/(N-1) );
y = 1./(alpha^2+x.^2);
[P,S] = polyfit(x,y,N-1);
yy_tcheby = polyval(P,xx);

plot( xx,yy,'k', xx,yy_lagrange,'k:', xx,yy_tcheby,'k-.' );
axis([-1,1,0,1/alpha^2+0.5]);
legend('Fonction à interpoler', 'Lagrange', 'Chebyshev');


saveas(gcf, '../../images/interpolation-chebyshev', 'eps');
saveas(gcf, '../../images/interpolation-chebyshev', 'png');
