%%
% Compare Taylor expansion with l2 expansion

addpath('../toolbox/');

name = 'l2';
name = 'taylor';
rep = MkResRep(name);

n = 1024;
t = linspace(-2*pi,2*pi,n)';

q = 25*2; r = ceil(50/q); it = 1;
for k=1:q
    s = (k-1)/(q-1);
    f = sin(t);
    
    [K,T] = meshgrid(0:k-1,t);
    D = T.^K;
    switch name
        case 'taylor'
            c = 1 ./ factorial(0:k-1)';
            c(1:2:end) = 0;
            c(2:2:end) = (-1).^(0:length(c(2:2:end))-1);
        case 'l2'
            c = pinv(D)*f;
    end
    
    clf; hold on;
    plot(t,f, 'k', 'LineWidth', 2);
    plot(t,D*c, 'LineWidth', 2, 'Color', [s 0 1-s]);
    axis([min(t) max(t) -2 2]);
    drawnow;
    set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1]); box on;
    for m=1:r
        % saveas(gcf, [rep 'anim-' znum2str(it,2) '.png' ]);
        it = it+1;
    end
end