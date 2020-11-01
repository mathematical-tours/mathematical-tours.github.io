function y = decode_hamming(x)
n = length(x); k = log2(n+1);
H = zeros(k,n); y = x;
for(j=1:n) H(:,j) = ecriture_binaire(j,k); end;
s = mod(H*x,2);
e = dot(s,2.^(0:(k-1)));
if(e~=0) y(e)=1-y(e); end;