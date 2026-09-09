clear;
clf;

N = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcule de la transformée de fourier d'un signal 
% qui n'est pas périodique
beta = N*0.23;
f = exp(2i*pi*(0:N-1)*beta/N);
f_prec = f;
f_prec(10*N) = 0;
ff = fft(f);
ff_prec = fft(f_prec);
xx_prec = 0:N/(10*N-1):N;
xx = 0:N-1;
plot( xx, abs(ff), 'k*', xx_prec, abs(ff_prec), 'k' );
legend('Spectre échantilloné', 'Spectre continu')
xlabel('Fréquence');
ylabel('Amplitude');
h = line( [beta beta], [0 20] );	% un trait au niveau du pic
set( h, 'LineStyle', ':' );
set(h,'Color',[0,0,0]);

saveas(gcf, '../spectre-signal-non-periodique', 'eps')
saveas(gcf, '../spectre-signal-non-periodique', 'png')

pause;
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% une fonction non périodique :
N = 128;
bruit = 0.2;
beta = 2.23;
f = cos(2*pi*(0:N-1)*beta/N) + bruit*(rand(1,N)-0.5);

% la fonction
subplot(2,2,1);
plot(0:N-1, f, 'k');
title('(a) Fonction analysée');
axis([0,N-1,-1.05,1.05]);
xlabel('Temps');
ylabel('Amplitude');

% sa transformée
ff = fft(f);
subplot(2,2,2);
plot(0:N/4-1, abs(ff(1:N/4)), 'k');	% le début du spectre
axis([0,N/4-1,0,60]);
title('(b) Transformée de Fourier');
xlabel('Fréquence');
ylabel('Amplitude');
h = line( [2 2], [0 100] );	% un trait au niveau du pic
set( h, 'LineStyle', ':' );
set(h,'Color',[0,0,0]);	

% zoom près du pic
delta = 1/N;
b = 2;
y = f .* exp( -2i*pi*b*(0:N-1)/N );
x = b:delta:(b+1-delta);
f_zoom = frft(y, delta);
subplot(2,2,3);
plot(x, abs(f_zoom), 'k');
axis([2,3,0,70]);
xlabel('Fréquence');
ylabel('Amplitude');
title('(c) Zoom sur le spectre');
h = line( [2.23 2.23], [0 100] );	% un trait au niveau du pic
set( h, 'LineStyle', ':' );
set(h,'Color',[0,0,0]);

% transformée réajustée
alpha = b/beta;
r = round(N*alpha);
ff_ajust = frft( f(1:r), r/(N*alpha) );
subplot(2,2,4);
plot(0:N/4-1, abs(ff_ajust(1:N/4)), 'k');	% le début du spectre
axis([0,N/4-1,0,60]);
title('(d) Transformée de Fourier ajustée');
xlabel('Fréquence');
ylabel('Amplitude');

saveas(gcf, '../transformee-fourier-fractionnaire', 'eps')
saveas(gcf, '../transformee-fourier-fractionnaire', 'png')