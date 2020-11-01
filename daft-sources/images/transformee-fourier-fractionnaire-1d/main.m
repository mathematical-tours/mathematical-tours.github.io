% exemples de transformée de Fourier fractionnaire (1D).

n = 128;
p = n/8;
x = [ones(p,1);zeros(n-p,1)];

subplot( 3,3,1 );
plot(0:n-1, x, '.-');
axis([0,n,-0.05,1.05]);
title('Fonction d''origine')

k=1;

for beta = 0.6:0.2:2.0

    k = k+1;
    
    y = frft( x, beta );
   % y = fftshift(y);
    subplot( 3,3,k );
    plot(0:n-1, real(y), '.-');
    axis tight;
    
    title(sprintf('\\alpha=%.2f', beta));

end

saveas(gcf, '../transformee-fourier-fractionnaire-1d', 'eps')
saveas(gcf, '../transformee-fourier-fractionnaire-1d', 'png')