
addpath('../toolbox/');
rep = MkResRep();


q = 100;
m = 10000; % #iter
m1 = 30; % displayed final point





mymap = @(x,r)r * ( sin(1.2*pi*x).^2 );
r_range = [0.2 1];


mymap = @(X,r)4*r*(X .* (1-X));
r_range = [.6 .99];

mymap = @(X,r)2*r*(.5 - abs(X-.5));
r_range = [.1 1];


mymap = @(x,r)4*r*(x .* (1-x)) + .1*sin(2*pi*x).^2;
r_range = [1 .99];


x = linspace(0,1,1000);
clf; plot(x, mymap(x,r_range(2)));

rlist = linspace(r_range(1),r_range(2), q); 




% precompute limit points 
clf; hold on;
for it=1:q    
    r = rlist(it);
    g = (it-1)/(q-1);
    
    x = 1/2;
    for i=1:m
        x(end+1) = mymap(x(end),r);
    end
    if 1
        s = ones(m1,1)*10; % size
        col = [g 0 1-g] .* ones(m1,1);
        scatter( r*ones(m1,1), x(end-m1+1:end)', s,  col, 'filled' );
        axis([r_range(1) 1.005*r_range(2) 0 1]);
        set(gca, 'XTick', 0:1:10);
    else
        nmax = 50;
        clf; hold on;
        plot(x(1:nmax), '.-', 'color', [g 0 1-g], 'LineWidth', 2, 'MarkerSize', 25);
        axis([0  nmax 0 1]);
        set(gca, 'XTick', 0:10:nmax);
    end
    set(gca, 'YTick', [0 .5 1], 'FontSize', 15);
    set(gca, 'PlotBoxAspectRatio', [1 2/3 1]);
	box on;
    drawnow;
    % saveas(gcf, [rep 'anim-' znum2str(it,3), '.png']);    
end
% drawnow;
    % imwrite(rescale(1-Xh), [rep 'anim-' znum2str(it,3) '.png']);
 
