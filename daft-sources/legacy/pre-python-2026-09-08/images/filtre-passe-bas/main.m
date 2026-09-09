% dessine la réponse fréquentielle du filtre (question 2)
N = 32; P = 1024;

leps = [0,0.3,0.7];
nb = length(leps);
c = [(N/2+1:N),(1:N/2)];
cc = [(P/2+1:P),(1:P/2)];
x = (-N/2):(N/2-1);
xx = ((-P/2):(P/2-1))*pi/(P-1);

for i=1:nb
    eps=leps(i);
    f = filtre_parametrable(N,eps);
    ff1 = real(fft(f));     % TFD
    ff2 = [f(1:N/2); zeros(P-N,1); f((N/2+1):N)];
    ff2 = real(fft(ff2));   % TF continue
    
    subplot(nb,3, (i-1)*3+1);
    plot(x,f(c),'k.:');
    axis tight;
    if( i==1 ) title('Filtre f'); end;
    ylabel( ['\epsilon = ', sprintf('%.1f',eps)] );
    
    subplot(nb,3, (i-1)*3+2);
    plot(x,ff1(c),'k.:');
    axis tight;
    if( i==1 ) title('TFD de f'); end;
    
    subplot(nb,3, (i-1)*3+3);
    plot(xx,ff2(cc),'k');
    axis tight;
    if( i==1 ) title('TF continue de f'); end;
end

saveas(gcf, '../../images/filtre-passe-bas', 'eps')
saveas(gcf, '../../images/filtre-passe-bas', 'png')