% Transforme un vecteur en nombre.
function res = vector2number(v,b)
N = length(v);
res = sum( v.*( b.^(0:N-1)' ) );