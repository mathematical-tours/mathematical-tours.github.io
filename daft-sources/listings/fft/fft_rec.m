% Algorithme FFT version récursive
function res = fft_rec(f,dir)
% N doit être une puissance de 2
N = length(f);
res = zeros(size(f));
if N==1    % fin de l'algorithme
    res(1) = f(1);
    return;
end
N1 = N/2;
% construction des deux sous-vecteurs
f_p = f(1:2:N);   
f_i = f(2:2:N);
% appels récursifs
f_p = fft_rec(f_p,dir);
f_i = fft_rec(f_i,dir);
% application de l'opérateur S
f_i = operateur_s(f_i,dir*0.5);
% mixage des deux résultats
res(1:N1)     = f_p + f_i;
res((N1+1):N) = f_p - f_i;