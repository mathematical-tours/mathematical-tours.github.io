function g = Approx(f,k)

n = size(f,1);
m = ceil(n/k);
g = f;
for i=1:k
    for j=1:k
        seli = (i-1)*m+1:min(i*m,n);
        selj = (j-1)*m+1:min(j*m,n);
        g(seli,selj) = mean(mean(f(seli,selj)));
    end
end


end