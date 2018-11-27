%%
% Plot of "Generalized extreme value distribution"

addpath('../toolbox/');
rep = MkResRep();

n = 2048;
A = 4;
s = linspace(-A,A,n);

q = 9;
xilist = linspace(-.7,.9,q);
%xilist = [-.5 0 .5];
q = length(xilist);

clf; hold on;
for i=1:q
    xi = xilist(i);
    t = (i-1)/(q-1);
    if not(xi==0)
        f = exp(-(1+xi*s).^(-1/xi)) .* ((xi*s+1)>0);
    else
        f = exp(-exp(-s));
    end
    h = diff(f); h = (n-1)/(2*A)*h; % /sum(h);
    plot(s(1:n-1),h, 'LineWidth', 2, 'Color', [t 0 1-t]);
end
%axis tight;
set(gca, 'XTick', -A:2:A, 'YTick', 0:.2:1, 'PlotBoxAspectRatio', [1 1/2 1], 'FontSize', 20);
axis([-A A 0 .6]); box on;
saveas(gcf, [rep 'extreme.eps'], 'epsc');