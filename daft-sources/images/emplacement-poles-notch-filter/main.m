N = 50;

poles = [0.9*exp( complex(0,pi/4) ), 0.9*exp( complex(0,-pi/4) )];
zzeros = [1.1*exp( complex(0,pi/4) ), 1.1*exp( complex(0,-pi/4) )];
x = (0:N)/N;
w = exp( 2i*pi/N );
cont = w.^(-(0:N));

% dessin des contours et des zéros.
plot( real(poles),imag(poles), 'k*', real(zzeros),imag(zzeros), 'k+', real(cont),imag(cont),'k' );
TITLE('Contours de calcul');
LEGEND('Poles', 'Zéros');
XLABEL('Re(z)');
YLABEL('Im(z)');
AXIS([-1.1,1.1,-1.1, 1.1]);

% sauvegarde l'image
saveas(gcf, '../emplacement-poles-notch-filter', 'eps')
saveas(gcf, '../emplacement-poles-notch-filter', 'png')