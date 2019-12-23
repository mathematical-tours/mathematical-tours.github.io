%%
% Just display Bernstein approximation polynomials


addpath('../toolbox/');
rep = MkResRep();

f = @(x)4*abs(x-.3) - 2*abs(x-.65);
f = @(x)abs(x-1/2);

n = 5;
x = linspace(0,1,1024)';

binom = @(n)factorial(n) ./( factorial(0:n) .* factorial(n:-1:0) );

% Just display basis polynomials
n = 6;
A = binom(n) .* x.^(0:n) .* (1-x).^(n:-1:0);
clf
plot(x, A); 
axis tight;


BernPoly = @(f,n,x)sum( f( (0:n)/n) .*  binom(n) .* x.^(0:n) .* (1-x).^(n:-1:0), 2);

nmax = 50;
nlist = [1:3 5 10 30];

clf; hold on;
plot(x,f(x), 'k', 'LineWidth', 2);
for it=1:length(nlist)
    n = nlist(it);
    s = (it-1)/(length(nlist)-1);
    plot(x,BernPoly(f,n,x), 'r', 'Color', [s 0 1-s]);
end
axis tight;
