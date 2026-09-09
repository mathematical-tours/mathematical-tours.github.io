% implémentation naive de la DCT.
function y = dct2_triv(x)
n = length(x);
y = zeros(n,1);
for j=0:n-1
    for k=0:n-1
        y(j+1) = y(j+1) + x(k+1)*cos( (k+1/2)*j*pi/n ); 
    end
end
