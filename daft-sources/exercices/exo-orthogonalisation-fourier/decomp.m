function y = decomp(f)

global g;

n = length(f);
y = zeros(n,1);

for k=0:n-1
    y(k+1) = dot( f,shift(g,k) );    
end