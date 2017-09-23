%%
% Plot the periodized sinc

R = 1024;
n = 16;
a = n;
x = linspace(-a,a,R)+1e-10;

% periodized
h  = 1/n * sin(pi*x) ./ sin(pi/n*x);

clf; hold on;
plot(x,sinc(x), 'r', 'LineWidth', 2);
plot(x,h, 'b', 'LineWidth', 2);
plot(-a:a,zeros(2*a+1,1), 'k.', 'MarkerSize', 25);
axis tight;
box on;