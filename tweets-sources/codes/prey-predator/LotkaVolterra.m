%%
% Prey-predator using different scheme 

n = 256*2; % for display

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

xlist = linspace(0.01,4,n);
ylist = linspace(0.01,2.5,n);
[Y,X] = meshgrid(ylist,xlist);
Z = X + 1i*Y;

a = 2/3; b = 4/3; gamma = 1; d = 1; 
xc = gamma/d; yc = a/b;  % center

f = @(x,y)a*x-b*x.*y;
g = @(x,y)d*x.*y-gamma*y;

V = @(x,y)d*x-gamma*log(x)+b*y-a*log(y);

K = 5;  % #particles
r = 20; % #levellines
m = linspace(0,1,r-1)';
CM = m*[1 .3 .3] + (1-m)*[.3 .3 1];
CM = m*[1 1 1] + (1-m)*[1 1 1]*.1;
col = distinguishable_colors(K);

% seed points
z0 = xc + 1i*yc + linspace(.1,1,K)'*.8*(1+1i);
x0 = real(z0); y0 = real(z0); 

tau = .05; 
tmax = 12;

meth = 'explicit';
meth = 'leapfrog';

niter = tmax/tau;
x = x0; y = y0; 
for j=1:niter
    switch meth
        case 'leapfrog'
            % leapfrog integrator
            x(:,end+1) = x(:,end) + tau/2 * f(x(:,end),y(:,end));
            y(:,end+1) = y(:,end) + tau * g(x(:,end),y(:,end));
            x(:,end) = x(:,end) + tau/2 * f(x(:,end),y(:,end));
        case 'explicit'
            % leapfrog integrator
            [x(:,end+1),y(:,end+1)] = deal( ...
                x(:,end) + tau * f(x(:,end),y(:,end)),  ...
                y(:,end) + tau * g(x(:,end),y(:,end)) ); 
    end
end


q = 100; % #display
disp_list = round(linspace(1,niter,q));
for it=1:q
    %
    j = disp_list(it);
    clf; hold on;
    imagesc(xlist,ylist,rescale(V(X,Y))');
    % contour(xlist,ylist,rescale(V(X,Y))',linspace(0,1,r), 'k');
    for k=1:K
        plot(x(k,1:j),y(k,1:j), 'color', col(k,:), 'LineWidth', 2);
        plot(x0(k), y0(k), '.', 'color', col(k,:), 'MarkerSize', 15);
    end
    plot(xc, yc, 'w.', 'MarkerSize', 15);
    colormap(CM);
    caxis([0 1]);
    axis equal; 
    axis image; axis off;
    % axis([0 max(xlist) 0 max(ylist)]);
    drawnow;
    mysaveas(it);
end


