clear;
N = 512;

x = (0:N-1)'*(0:N-1);
% matrices de Fourier
F = 1/sqrt(N)*exp(2*i*pi/N).^x;
% matrice intermédiaire pour le calcul des vecteurs propres
x = 0:N-1;
d = 2*( cos(2*pi/N*x) - 2 );
S = diag(d,0) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
S(1,N) = 1;
S(N,1) = 1;

[V,D] = eig(S);

DF = V^(-1)*F*V;

AF = mod(ceil(angle(diag(DF))*2/pi), 4);

A0 = [];
A1 = [];
A2 = [];
A3 = [];
% remise dans l'ordre
for i=1:N
    if AF(i)==0
        A0 = [A0, V(:,i)];
    elseif AF(i)==1
        A1 = [A1, V(:,i)];
    elseif AF(i)==2
        A2 = [A2, V(:,i)];
    else
        A3 = [A3, V(:,i)];
    end
end

A = [A0,A1,A2,A3];

imagesc(-abs(V));
colormap gray;
axis image;
set(gca, 'XTick', []);
set(gca, 'YTick', []);

% title('Matrice de vecteurs propres orthogonaux de la TFD');


saveas(gcf, '../../images/matrice-vecteurs-propres-tfd', 'eps');
saveas(gcf, '../../images/matrice-vecteurs-propres-tfd', 'png');

% subplot(2,2,1);
subx = 3*N/8:(5*N/8);
plot( subx, abs(V(subx,1)), 'k', subx, abs(V(subx,2)), 'k:', subx, abs(V(subx,3)), 'k--'  );
m = max(abs(V(subx,1)));
axis([3*N/8,(5*N/8),0,1.05*m]);
% title('Les trois premiers vecteurs propres');
legend('premier', 'deuxième', 'troisième');


saveas(gcf, '../../images/vecteurs-propres-tfd', 'eps');
saveas(gcf, '../../images/vecteurs-propres-tfd', 'png');