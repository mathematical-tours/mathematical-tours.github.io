% just generate signals and mixing

addpath('../toolbox/');
rep = MkResRep();

myset = @()set(gca, 'PlotBoxAspectRatio', [1 1/2 1], 'FontSize', 20, 'XTick', [], 'YTick', [], 'box', 'on');
mysaveas = @(name)saveas(gcf, [rep name '.eps'], 'epsc');

n = 256*2;
t = linspace(0,1,n);

a = sin(2*pi*7*t);
b = 2*mod(t*6.25,1)-1; 

plot(a,b,'.');


v = 1*exp(1i*pi*.25);
w = .6*exp(1i*pi*.40);
A = a*real(v) + b*real(w);
B = a*imag(v) + b*imag(w);






clf; plot(t,a, 'r', 'LineWidth', 2); myset(); axis tight;
mysaveas('a-init');
clf; plot(t,b, 'b', 'LineWidth', 2); myset(); axis tight;
mysaveas('b-init');

clf; plot(t,A, 'color', [0 1 0]*.6, 'LineWidth', 2); myset(); axis tight;
mysaveas('a-mixed');
clf; plot(t,B, 'color', [1 0 1]/2, 'LineWidth', 2); myset(); axis tight;
mysaveas('b-mixed');

clf; plot(a, b, 'k.', 'MarkerSize', 15); axis equal; axis off;
mysaveas('points-init');
clf; plot(A, B, 'k.', 'MarkerSize', 15); axis equal; axis off;
mysaveas('points-mixed');