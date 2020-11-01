function y = haar(x)
n = length(x); j = log2(n); y = [];
for i=0:j-1
    d = 1/sqrt(2)*( x(1:2:2^(j-i))-x(2:2:2^(j-i)) ); 
    y = [y;d];
    x = 1/sqrt(2)*( x(1:2:2^(j-i))+x(2:2:2^(j-i)) );
end
y = [y;x];