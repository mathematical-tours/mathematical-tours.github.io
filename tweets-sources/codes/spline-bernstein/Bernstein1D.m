%%
% Just display Bernstein approximation polynomials


addpath('../toolbox/');
rep = MkResRep();

f = @(x)4*abs(x-.3) - 2*abs(x-.65);
f = @(x)abs(x-1/2)<.2;

f = @(x)abs(x-1/2);
f = @(x)sin(7*pi*x) - 6*(x-1/2).^2;

m = 1024;
x = linspace(0,1,m)';

binom = @(n)factorial(n) ./( factorial(0:n) .* factorial(n:-1:0) );

% Just display basis polynomials
n = 6;
A = binom(n) .* x.^(0:n) .* (1-x).^(n:-1:0);
clf
plot(x, A); 
axis tight;


BernPoly = @(f,n,x)sum( f( (0:n)/n) .*  binom(n) .* x.^(0:n) .* (1-x).^(n:-1:0), 2);

q=70;
nlist = unique(1+round(linspace(0,1,q).^7 * 300));
%
clf; hold on;
plot(x,f(x), 'k', 'LineWidth', 2);
for it=1:length(nlist)
    n = nlist(it);
    s = (it-1)/(length(nlist)-1);
    % y = BernPoly(f,n,x);
    y = BernsteinBasis(x,n+1) * f( (0:n)'/n );
    %
    plot(x,y, 'r', 'Color', [s 0 1-s]);
    axis([0 1 min(f(x))-.05 max(f(x))+.05]); 
    box on;set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1]);
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
    drawnow;
end

for it=1:length(nlist)
    n = nlist(it);
    s = (it-1)/(length(nlist)-1);
    %
    clf;
    plot( x, BernsteinBasis(x,n+1) );
    axis([0 1 0 1]); 
    box on;set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1]);
    drawnow;
    saveas(gcf, [rep 'basis-' znum2str(it,2) '.png']);
end