n = 16; nn = 200; alpha = 0.3;
x = cos( ((0:n-1)+1/2)*pi/n )';
f = 1./(alpha^2+x.^2);
coef = 2/n*dct2(f);
coef(1) = coef(1)*1/2;
xx = (-1:2/(nn-1):1)';
ff = zeros(nn,1);
for k=0:n-1
    ff = ff + coef(k+1)*cos(k*acos(xx));
end
plot(x,f,'o',xx,ff);