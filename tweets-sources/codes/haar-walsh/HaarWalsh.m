%%
% Contrast Haar and Walsh transforms.

rep = '../results/haar-walsh/';
[~,~] = mkdir(rep);

n = 256;
t = (0:n-1)'/n-.5;
x = cos(4*pi*t);

x = .5*(abs(2*t+.5)<.2) - (abs(2*t-.3)<.2)*.8;


SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);



j = log2(n);
Nlist = n ./ 2.^(0:j);
Nlist = n ./ [1 2 4 n];

for j=1:length(Nlist)
    nmax = Nlist(j);
    t = (j-1)/(length(Nlist)-1);
    %%
    y = Haar(x,nmax);
    clf;
    plot( y, 'Color', [t 0 1-t], 'LineWidth', 2 );
    box on; axis tight; 
    axis([1 n -max(abs(y)) max(abs(y))]); 
    SetAR(1/3); set(gca, 'XTick', [], 'YTick', []);
    saveas(gcf, [rep 'haar-' num2str(nmax) '.eps'], 'epsc');
    %%
    y = Walsh(x,Nlist(j));
    clf;
    plot( y, 'Color', [0 t 1-t], 'LineWidth', 2 );
    box on; axis tight;
    axis([1 n -max(abs(y)) max(abs(y))]); 
    SetAR(1/3); set(gca, 'XTick', [], 'YTick', []);
    saveas(gcf, [rep 'walsh-' num2str(nmax) '.eps'], 'epsc');    
end