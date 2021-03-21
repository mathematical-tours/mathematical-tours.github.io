% create base points

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

% square points
p = 8;
z = exp(2i*pi*(0:p-1)'/p);
z = .45*z ./ max( abs(real(z)),abs(imag(z)) );

q = 75;
theta_list = linspace(0,2*pi,q);

radius = .05; % turn radius
dt = .01/6;

for it=1:q
    
    theta = theta_list(it);
    
    s = (it-1)/(q-1);
    s = (cos(2*pi*s)+1)/2;
    col = [1-s 0 s];
    
    clf; hold on;
    for k=1:p
        p1 = [0,0,0];
        p2 = [real(z(k)),imag(z(k)),theta];
        g = dubins_curve(p1,p2, radius, dt, true);
        plot(g(:,1), g(:,2), 'color', col, 'LineWidth', 2); 
        sc = .2;
        plot([p1(1),p1(1) + sc*cos(p1(3))],[p1(2),p1(2) + sc*sin(p1(3))], 'color', 'k', 'LineWidth', 2);
        plot([p2(1),p2(1) - sc*cos(p2(3))],[p2(2),p2(2) - sc*sin(p2(3))], 'color', col, 'LineWidth', 2);
        plot(p1(1),p1(2), '.', 'color', 'k', 'MarkerSize', 20);
        plot(p2(1),p2(2), '.', 'color', col, 'MarkerSize', 20);
    end
    axis equal; axis([-1 1 -1 1]); box on;
    set(gca, 'XTick', [], 'YTick', [])
    drawnow;
    mysaveas(it);
end