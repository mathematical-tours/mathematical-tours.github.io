% nombre de points d'interpolation pour le calcul de l'intégrale
M = 2^8; h = 1/M;
% valeur de f aux points d'interpolation
f_val = zeros(M,1);
for i=1:M
    f_val(i) = f((i-1)*h);
end
% calcul de la fft
dft_val = fft(f_val);
% calcul des coefficients de fourier : 
%  - renverser le tableau (à cause du '-' dans le -2*i*pi des coef de Fourier)
%  - multiplier par h
fcoef = zeros(M,1);
for n=1:M
    i = 1 + mod(-(n-1),M);
    fcoef(n) = h*dft_val(i);
end
% dessine une évolution de la solution
t = [0.01,0.02,0.03,0.04,0.05,0.06]';
ts = length(t);
xx = real( solve_eq(0, fcoef) );
for i=1:ts
    xx = [xx, real( solve_eq(t(i), fcoef) )];
end
plot(xx);