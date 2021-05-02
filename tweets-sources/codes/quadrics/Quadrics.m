addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name, it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);

xmax = 3;
n = 40;
x = linspace(-xmax,xmax,n);
[X,Y,Z] = ndgrid(x,x,x);

a = 1;
b = -1;

q = 60*2;
for it=1:q
    t = (it-1)/(q-1);
    F = X.^2 + cos(2*pi*t)*Y.^2 + sin(2*pi*t)*Z.^2;
    T = 1; % level
    col = [.5+.5*cos(2*pi*t) 0 .5+.5*sin(2*pi*t)];
    clf;
    p = patch( isosurface( x,x,x,F,T ) );
    isonormals( x,x,x,F,p );
    set(p, 'FaceColor', col, 'EdgeColor', 'none');
    % isocolors(x,x,x,R(:,:,:,1),R(:,:,:,2),R(:,:,:,3),p);
    % p.FaceColor = 'interp';
    % p.EdgeColor = 'none';
    box on; axis on;
    set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    axis equal;
    axis([-1 1 -1 1 -1 1]*xmax);
    lighting gouraud;
    view(3); axis off;
    camlight; 
    drawnow;
    mysaveas('anim', it);
end


q = 60*2;
m = 1.3;
for it=1:q
    t = (it-1)/(q-1);
    col = [.5+.5*cos(2*pi*t) 0 .5+.5*sin(2*pi*t)];
    clf; hold on;
    plot([-m m], [0 0], 'k:', 'LineWidth', 1);
    plot([0 0], [-m m], 'k:', 'LineWidth', 1);
    axis equal; axis([-1 1 -1 1]*m); 
    box on;
    set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    plot(cos(2*pi*t), sin(2*pi*t), '.', 'MarkerSize', 45, 'col', col);
    drawnow;
    mysaveas('dot', it);
end
