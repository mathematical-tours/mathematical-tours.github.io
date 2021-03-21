%%
% Polynomial interpolation

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


f = @(x)1./(1+25*x.^2);

t = linspace(-1,1,500)';



it = 1;
imax = 30;
for i=1:imax
    n = i;
    s = (i-1)/(imax-1);
    %
    xi = linspace(-1,1,n)';
    xi = cos( pi*linspace(0,1,n) );
    yi = f(xi);
    %
    c = pinv( xi(:) .^ (0:n-1) ) * yi(:);
    g = ( t(:) .^ (0:n-1) ) * c;
    %
    clf; hold on;
    plot(t,f(t), 'k', 'LineWidth', 2);
    plot(t,g, 'color', [s 0 1-s], 'LineWidth', 2);
    plot(xi,f(xi), '.', 'color', [s 0 1-s], 'MarkerSize', 25);
    box on;
    axis([-1 1 -.5 1.5]);
    set(gca, 'XTick', [], 'YTick', []); 
    drawnow;
    for j=1:4
        mysaveas(it);
        it = it+1;
    end
end