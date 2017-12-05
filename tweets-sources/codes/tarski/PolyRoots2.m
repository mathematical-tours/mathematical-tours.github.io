%%
% plot set of polynomials with real roots

n = 200;
a = linspace(-1,1,n);
b = linspace(-1,1,n);
x = linspace(-1,1,n);

T = zeros(n,n);
for i=1:n
    for j=1:n
        r = roots([1, 0, a(i), 0, b(j)]);
        T(i,j) = sum(abs(imag(r))<1e-4)>0;
    end
end

n = 120;
a = linspace(-1,1,n);
b = linspace(-1,1,n);
x = linspace(-1,1,n);

[A,B,X] = ndgrid(a,b,x);
U = X.^4 + A.*X.^2 + B;

options.isolevel = 0;
options.alpha = 1;
options.color = 'blue';
clf;
F = plot_isosurface(U,options);
axis on; box on;
view(145,50)
