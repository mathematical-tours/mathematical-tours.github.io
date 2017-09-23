%%
% plot of Dirichlet kernel -> Dirac

Nlist = [1 4 6 10 15];

x=linspace(-2*pi,2*pi,1024)'+1e-9;

D = [];
clf; hold on;
for i=1:length(Nlist)
    n = Nlist(i);
    D = 1/(n+1/2) * sin( (n+1/2)*x ) ./ sin(x/2);
    t = (i-1)/(length(Nlist) - 1);
    plot(x,D,'Color', [1-t 0 t], 'LineWidth', 2);
end
axis tight; box on;