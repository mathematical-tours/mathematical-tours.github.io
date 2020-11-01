% Transforme un nombre en vecteur.
function res = number2vector(x,b)
N = floor( log(x)/log(b) )+1;
res = zeros(N,1);
for i=1:N
	q = floor(x/b);
	res(i) = x - q*b;
	x = q;
end