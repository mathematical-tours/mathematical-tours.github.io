% Opérateur Chi
function res = operateur_chi(a,x)
N = length(a);
a_inv = [a(1);a(N:-1:2)];
res = a.*cos( 2*pi*x*(0:N-1)'/N ) + a_inv.*sin( 2*pi*x*(0:N-1)'/N );