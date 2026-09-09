% précision du tracé
N = 50;
h = 1/N;
zt = zeros(N+1,N+1);

for i=1:N+1
    x = (h*i)*3 - 1.5;
    for j=1:N+1
        y = (h*j)*3 - 1.5;
        z = complex(x,y);
        zt( i,j ) = fonction_trans(z);
    end
end


xx = (0:h:1)*3 - 1.5;
clf;
surf(xx, xx, abs(zt));
AXIS([-1.5,1.5,-1.5,1.5,0.5,1.5]);
VIEW(120,40);
CAXIS([0.5,1.5]);
XLABEL('X');
YLABEL('Y');
ZLABEL('Amplitude');
saveas(gcf, '../notch-filter-trans-Z', 'eps')
saveas(gcf, '../notch-filter-trans-Z', 'png')

pause;

% répsonse fréquencielle
N = 256;
h = 1/N;
rf = zeros(N,1);

for i=1:N
    t = (i-1)*h*2*pi;
    z = exp( complex(0,t) );
    rf( i ) = fonction_trans(z);
end
xx = 0:h:(1-h);
subplot(2,1,1);
plot(xx, abs(rf), 'k');
XLABEL('Fréquences');
YLABEL('Amplitude');
AXIS([0,0.5,0, 1.2]);
TITLE('Réponse fréquencielle')

% réponse impulsionnelle 
% (la réponse fréquencielle est symétrique, la fonction de transfert est réelle)
ri = real( ifft(rf) );
subplot(2,1,2);
plot( 0:N-1, ri, 'k:*' );
AXIS([0,20,-0.2, 1]);
TITLE('Réponse impulsionnelle');
XLABEL('Temps');
YLABEL('Amplitude');

% sauvegarde l'image
saveas(gcf, '../notch-filter-rep-imp-freq', 'eps')
saveas(gcf, '../notch-filter-rep-imp-freq', 'png')