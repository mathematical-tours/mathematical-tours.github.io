function g = ApproxInterp(f,k)


n = size(f,1);
I = round(linspace(0,n,k+1));

for i=1:k
    for j=1:k
        seli = I(i)+1:I(i+1);
        selj = I(j)+1:I(j+1);
        for s=1:size(f,3)
            g(i,j,s) = mean(mean(f(seli,selj,s)));
        end
    end
end

% re-sample

U = floor( linspace(1,k+1,n+1) ); U = U(1:n);
g = g(U,U);

end

