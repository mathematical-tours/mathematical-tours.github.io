
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


radius = .25; % turn radius
dt = .01;
L = @(q)sum( sqrt( sum( (q(2:end,1:2)-q(1:end-1,1:2)).^2, 2) ) );



n = 101;
x = linspace(-1,1,n);
D = zeros(n);

q = 75;
theta_list = linspace(0,2*pi,q);
r = 15; % #levellines

for it=1:q
    theta = theta_list(it);
    for i=1:n
        for j=1:n
            D(i,j) = L( dubins_curve([0,0,0],[x(i),x(j),theta], radius, dt, true) );
        end
    end
    s = (it-1)/(q-1);
    s = (cos(2*pi*s)+1)/2;
    m = linspace(0,1,r-1)';
    CM = m*[1-s 0 s] + (1-m)*[1 1 1];
    % display
    clf; hold on;
    imagesc(x,x,D);
    contour(x,x,D,linspace(0,2,r), 'k');
    colormap(CM);
    caxis([0 2]);
    axis image; axis off;
    drawnow;
    mysaveas(it);
end




