CLF;
% Montre la différence entre la convolution 
% linéaire et la convolution linéaire.
N = 16;
f = [0:N/2, (N/2-1):-1:1]';
P = 3;

g = [ones(P,1);zeros(N-2*P,1);ones(P,1)];
gg = [ones(P,1);zeros(2*N-2*P,1);ones(P,1)];
sel = [((N/2+1):N),1:N/2];

subplot(2,2,1);
plot(0:N-1, f, 'k.:');
title('f');

subplot(2,2,2);
plot((-N/2+1):N/2, g(sel), 'k.:');
axis tight;
title('g');

a = real( ifft( fft(f).*fft(g) ) );

subplot(2,2,3);
plot(0:N-1, a, 'k.:');
title('Convolution cyclique');

ff = [f;zeros(N,1)];
b = real( ifft( fft(ff).*fft(gg) ) );
b = b(1:N);

subplot(2,2,4);
plot(0:N-1, b, 'k.:');
title('Convolution acyclique');

saveas(gcf, '../diff-conv-lineaire-circulaire', 'eps')
saveas(gcf, '../diff-conv-lineaire-circulaire', 'png')