% Opérateur S
function res = operateur_s(a,x)
N = length(a);
res = a.*exp( -2.0i*x*(0:N-1)'*pi/N );