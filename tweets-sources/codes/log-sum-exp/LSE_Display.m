%%
% Check for smoothed convex enveloppes

addpath('../toolbox/');
rep = MkResRep();


n = 1024; % for display
x = linspace(-1,1,n)';

% sample points
m = 10;

xi = linspace(-.8,.8,m);
yi = 1/2 * xi.^2 + .15*xi;
yDi = xi + .15;



M = yi + (x-xi) .* yDi;

smax = @(M,epsilon)epsilon*log(sum(exp(M/epsilon),2));



q = 50;

elist = linspace(.6,.001,q/2);
elist = linspace(.3,.001,q/2);

elist = elist(floor(1:.5:q/2+.5));

for i=1:q
    clf; hold on;
    for j=1:i
        s = (j-1)/(q-1);
        epsilon = elist(j);
        col = [s 0 1-s];
        if i~=j
            col = .2*col + .8;
        end
        plot(x, smax(M,epsilon), 'color', col, 'LineWidth', 2);
    end
    plot(x, M, 'color', [1 1 1]*.3);
    axis([-1 1 -.3 .6]);
    box on;
    set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'XTick', [], 'YTick', []);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(i,2) '.png']);
end






