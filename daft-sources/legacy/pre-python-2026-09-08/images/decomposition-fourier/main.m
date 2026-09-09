% premier contact avec la décomposition de Fourier.
N = 17;
N0 = 8;

c = [((N0+2):N)';(1:N0+1)'];

ramp	= [-N0:N0];
ramp	= pi * ramp/N;
f	= exp( -(ramp.^2)/(0.5) );
f = f(c);

ff = real( fft(f) );


subplot(1,2,1);
plot(0:N-1,f(c),'k*:');
axis square;
axis tight;
title('Fonction originale');

subplot(1,2,2);
plot(-N0:N0,ff(c),'k*:');
axis square;
axis tight;
title('Transformée de Fourier');


saveas(gcf, '../decomposition-fourier', 'eps');
saveas(gcf, '../decomposition-fourier', 'png');