% Transformée de Hartley rapide
function res = fht(f)
% N doit être une puissance de 2
N = length(f); N1 = N/2;
res = zeros(size(f));
if N==1    % fin de l'algorithme
    res(1) = f(1);
    return;
end
% construction des deux sous-vecteurs
f_p = f(1:2:N); f_i = f(2:2:N);
% appels récursifs
f_p = fht(f_p); f_i = fht(f_i);
% application de l'opérateur chi
f_i = operateur_chi(f_i,0.5);
% mixage des deux résultats
res(1:N1)     = f_p + f_i;
res((N1+1):N) = f_p - f_i;