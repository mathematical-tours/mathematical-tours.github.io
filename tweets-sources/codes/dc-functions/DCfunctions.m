n = 1024;

t = (0:n-1)'/n;
R = @(x)x-mean(x(:));
%R = @(x)x;


rep = ['results/'];
[~,~] = mkdir(rep);

name = 'affine';
name = 'gaussian';
name = 'oscillations';

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);


switch name
    case 'oscillations'
        f = cos(2*pi*4*t) + 10*exp(-(t-.4).^2/.5);
    case 'gaussian'
        f = exp(-(t-.5).^2/.05);
    case 'affine'
        f = abs(t-.3) - .6*abs(t-.5) + .5*abs(t-.7);
end
f = R(f);


A = R(cumsum(R(cumsum(max(diff(f,2),0)))));
B = R(cumsum(R(cumsum(min(diff(f,2),0)))));

clf;
plot(f, 'k', 'LineWidth', 2); axis tight; 
set(gca, 'XTick', [], 'YTick', []);
SetAR(2/3);
saveas(gcf, [rep name '-func.eps'], 'epsc');

clf;
plot([A B], 'LineWidth', 2); axis tight; 
set(gca, 'XTick', [], 'YTick', []);
SetAR(2/3);
saveas(gcf, [rep name '-dc.eps'], 'epsc');

