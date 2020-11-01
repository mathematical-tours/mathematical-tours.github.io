% calcul "à la main" la réponse impulsionnelle.
a0 = 1; a1 = -1.414; a2 = 1;
b1 = 1.273; b2 = -0.810;

N = 200;
ri = zeros(N+2,1);
x = zeros(N+2,1); x(3) = 1; 
for i=3:(N+2)
	ri(i) = a0*x(i) + a1*x(i-1) + a2*x(i-2) + b1*ri(i-1) + b2*ri(i-2);
end

ri = ri(3:(N+2));
rf = real(fft(ri));


% calcul des contours
w = exp( 2i*pi/N );
a = 1;
cont1 = a * w.^(-(0:N));
y1 = czt(ri,N+1,w,a);

poles = [0.9*exp( complex(0,pi/4) ), 0.9*exp( complex(0,-pi/4) )];
zzeros = [1.1*exp( complex(0,pi/4) ), 1.1*exp( complex(0,-pi/4) )];

w = 1/(1 + 1/(10*N))*exp( 2i*pi/N );
a = 1;
cont2 = a * w.^(-(0:N));
y2 = czt(ri,N+1,w,a);

x = (0:N)/N;

% dessin des contours et des zéros.
subplot(1,2,1);
plot( real(poles),imag(poles), 'k*', real(zzeros),imag(zzeros), 'k+', real(cont1),imag(cont1), 'k', real(cont2),imag(cont2), 'k-.' );
title('Contours de calcul');
axis square;
legend('Poles', 'Zéros');
xlabel('Re(z)');
ylabel('Im(z)');
axis([-1.2,1.2,-1.2,1.2]);

% dessin des transformées
subplot(1,2,2);
plot( x, abs(y1), 'k' ,  x, abs(y2), 'k-.' );
TITLE('Transformée en Z sur le contour.');
XLABEL('Fréquence');
YLABEL('Amplitude');
legend('contour n°1', 'contour n°2');
axis square;

% sauvegarde l'image
saveas(gcf, '../z-transform', 'eps')
saveas(gcf, '../z-transform', 'png')