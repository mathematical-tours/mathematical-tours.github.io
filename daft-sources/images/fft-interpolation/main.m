N = 11;

x = 0:N-1;
f = rand(N,1);

[xx,ff] = fft_interp(f,100);

plot(x,f,'k',xx,ff,'k--');
legend('Courbe originale', 'Courbe interpolée');

saveas(gcf, '../fft-interpolation', 'eps');
saveas(gcf, '../fft-interpolation', 'png');

