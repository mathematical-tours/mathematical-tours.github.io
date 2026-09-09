CLF;
% nombre de points d'interpolation pour le calcul de l'intégrale [méthode des rectangles à 
% gauche, pour que la somme devienne une FFT]
M = 2^8;
h = 1/M;
% NOTE : si on était intéressé par des pulsations autres que 2*n*Pi (ie.les coef de fourier), 
% on pourrait utiliser une astuce, en formant un vecteur f_val de longueur N>M, avec des 0 à la
% fin, et écrire que l'on veut 2*pi*n/M sous la forme 2*pi*p*M/N [pour calculer le
% pième coef d'une FFT de longueur N].
% mais alors on a n=p*M/N qui est non-entier si M!=N (mais ceci permet de calculer des valeurs 'intermédiaires')

% valeur de f aux points d'interpolation
% vaut 0 pour les indices > M-1
f_val = zeros(M,1);
for i=1:M
	f_val(i) = f((i-1)*h);
end

% calcul de la fft
dft_val = zeros(M, 1);
dft_val = fft(f_val);

% calcul des coef de fourier : renverser le tableau (à cause du '-' dans le -2*i*pi des coef de Fourier)
% et multiplier par h
% attention, les coefs sont stockés dans l'ordre de la FFT : ie 0,1,...,M/2,-M/2+1,...,-1
fcoef = zeros(M,1);
for n=1:M
   i = 1 + mod(-(n-1),M);
	fcoef(n) = h*dft_val(i);
end

% dessine une évolution de la solution
t = [0.01,0.02,0.03,0.04,0.05,0.06]';
[ts,tss] = size(t);
%xx = zeros(301,1);
xx = real( solve_eq(0, fcoef) );

% for i=1:ts
% 	xx = [xx, real( solve_eq(t(i), fcoef) )];
% end
% plot(xx);

N = 301;
h = 1/300;
for i=0:300
    f_dep(i+1) = f(i*h);
end
s1 = real( solve_eq(0.001, fcoef) );
s2 = real( solve_eq(0.01, fcoef) );
s3 = real( solve_eq(0.1, fcoef) );
xx = 0:h:1; 
plot(   xx, f_dep, 'k:', xx, s1, 'k--', xx, s2, 'k-.', xx, s3, 'k-');
AXIS([0,1,0,1.1]);
XLABEL('X');
LEGEND('t=0', 't=0.001', 't=0.01', 't=0.1');

saveas(gcf, '../eq-chaleur-evolution', 'eps')
saveas(gcf, '../eq-chaleur-evolution', 'png')