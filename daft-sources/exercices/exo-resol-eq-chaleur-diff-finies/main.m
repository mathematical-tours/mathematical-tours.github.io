n = 100;
x = (0:1/(n-1):1)';
x = cos(2*pi*x) + rand(n,1)*0.2;
% x = resolution_explicite(x,0.2,4);
x = resolution_implicite(x,0.2,0.5,5);
plot(x);