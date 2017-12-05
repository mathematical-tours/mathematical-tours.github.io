%%
% plot set of polynomials with real roots

n = 60;
a = linspace(-1,1,n);
b = linspace(-1,1,n);
c = linspace(-1,1,n);


T = zeros(n,n,n);
for i=1:n
    for j=1:n
        for k=1:n
            r = roots([1, 0, a(i), b(j), c(k)]);
            T(i,j,k) = sum(abs(imag(r))<1e-4)>0;
        end
    end
end

options.isolevel = 1/2;
options.alpha = 1;
clf;
F = plot_isosurface(T,options);

n = 120;
a = linspace(-1,1,n);
b = linspace(-1,1,n);
c = linspace(-1,1,n);

[A,B,C] = ndgrid(a,b,c);
Delta = 16*A.^4 .* C - 4*A.^3 .* B.^2 - 128*A.^2 .* C.^2 + 144*A.*B.^2 .* C - 27*B.^4 + 256*C.^3;
options.isolevel = 0;
options.alpha = 1;
clf;
F = plot_isosurface(Delta,options);