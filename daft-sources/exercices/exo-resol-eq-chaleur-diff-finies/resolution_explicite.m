function y = resolution_explicite(x,s,t)
n = length(x);  y = x;
for i = 1:t
    for k=1:n
        y1(k) = s*y( mod(k-2,n)+1 ) + s*y( mod(k,n)+1 ) + (1-2*s)*y(k);
    end
    y = y1;
end
