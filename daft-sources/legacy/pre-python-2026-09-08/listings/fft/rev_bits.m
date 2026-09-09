% Classe les élements du vecteur.
function y = rev_bits(x)
n = length(x);
t = floor(log2(n));
y = zeros(n,1);
for i=0:n-1
    j = rev_index(t,i);
    y(j+1) = x(i+1);
end