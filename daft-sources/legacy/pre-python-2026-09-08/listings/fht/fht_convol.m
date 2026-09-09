% convolution par FHT
function res = fht_convol(x,y)
% N doit être une puissance de 2
N = length(x);
res = zeros(size(x));
a = fht(x);
b = fht(y);
a_inv = [a(1);a(N:-1:2)];
b_inv = [b(1);b(N:-1:2)];
res = 0.5*( a.*b - a_inv.*b_inv + a.*b_inv + a_inv.*b );
res = fht(res)/N;