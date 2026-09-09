% test de l'interpolation trigonométrique.

n = 33;
eta = 5;
p = eta*n;
x = rand(n,1);
y = interp_trigo(x,eta);

a = 0:(n-1);
b = 0:n/p:(n-1/p);
plot(a,x,'o', b,y);