function g = approx_constant(f,Q)

[N,P] = size(f);

n = round(N/Q);
R = round(P/n);
g = f*0;


for i=1:ceil(N/n)
    seli = (i-1)*n+1:i*n;
    if i==ceil(N/n)
        seli = (i-1)*n+1:N;
    end
    for j=1:ceil(P/n)
        selj = (j-1)*n+1:j*n;
        if j==ceil(P/n)
            selj = (j-1)*n+1:P;
        end
        g(seli,selj) = mean(mean(f(seli,selj)));
    end
end



end
