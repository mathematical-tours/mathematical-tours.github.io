function y = ecriture_binaire(x,k)
y = zeros(k,1);
for(i=1:k) q = floor(x/2); y(i) = x-2*q; x = q; end;