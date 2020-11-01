N = 10; P = 512;
x = rand(N,1);

xx = (0:P-1)'/(P-1)*(N-1);

c = coef_spline_1(x);

y = zeros(P,1);
for k=0:N-1
    y = y + c(k+1)*spline3(xx-k);
end

plot(0:N-1, x, 'k*', xx, y, 'k');