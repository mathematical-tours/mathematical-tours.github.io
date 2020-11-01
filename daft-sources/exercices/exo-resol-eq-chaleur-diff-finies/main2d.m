n = 100;
x = (0:1/(n-1):1)'*(0:1/(n-1):1);
x = cos(4*pi*x) + rand(n,n)*0.2;
x = resolution_implicite_2d (x,0.2,0.5,5);
image(255*abs(x));