%%
% test on Legendre transforms

% f^*(y) = max xy - f(x)

f = @(x)1/2*x.^2;
f = @(x).4*abs(x) + 1/2*x.^2;

y1 = linspace(-1,1,512);
x1 = linspace(-10,10,1024*100);
[Y,X] = meshgrid(y1,x1);


x = linspace(-1,1,1023);
clf; 
plot(x,f(x), 'b-', 'LineWidth', 2);


[fs,I] = max(X.*Y-f(X));
xs = x1(I);

clf; 
plot(y1,fs, 'r-', 'LineWidth', 2);