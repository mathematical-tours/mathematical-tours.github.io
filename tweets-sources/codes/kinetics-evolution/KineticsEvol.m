%%
% Display of momentum equation in phase space

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

xmax = 1.7;
vmax = 1.7;
f = @(x)1/2 * (x.^2-1).^2;
nablaf = @(x)2*x.*(x.^2-1);

t = linspace(-xmax,xmax,300);

clf; hold on;
plot(t,f(t), 'LineWidth', 2);
% plot(t,nablaf(t), 'r');
box on;
set(gca, 'XTick', [], 'YTick', []);
axis([-xmax xmax 0 1.5]);
saveas(gcf, [rep 'fonction.eps'], 'epsc');

% beta*x'' + (1-beta)*x' = -nablaf(x)
% x' = v
% v' = 1/beta * ( - (1-beta)*x' -nablaf(x) )
%    = (beta-1)/beta*v - 1/beta*nablaf(x) =: W(x,v) 



beta = 0.2;  
dt = .005;


beta = 0.8;  
dt = .02;






beta = 0.01; 
dt = .02/4;




beta = 1;  
dt = .015;

beta = .9;  
dt = .02;

beta = 0.5;  
dt = .01/2;

W = @(x,v)(beta-1)/beta*v - 1/beta*nablaf(x);

disp_mode = 'density';
disp_mode = 'particles';

switch  disp_mode
    case 'particles'
        x = rescale(rand(p,1),-m*xmax,m*xmax);
        v = rescale(rand(p,1),-m*vmax,m*vmax);
        p = 200;
    case 'density'
        p = 60000; % #particles
        p0 = round(sqrt(p)); p = p0^2;
        m = .5;
        [v,x] = meshgrid(linspace(-m*vmax,m*vmax,p0),linspace(-m*xmax,m*xmax,p0));
        x = x(:);
        v = v(:);
end



niter = 2000;

k = 1;
xvg = x;
for it=1:niter
    if mod(it,10)==1
        
        clf; hold on;
        if strcmp(disp_mode,'density')
            M = 1500;
            xlist = linspace(-xmax,xmax,M);
            vlist = linspace(-vmax,vmax,M);
            h = parzen2d( (x+xmax)/(2*xmax), (v+vmax)/(2*vmax),M,15);
            imagesc(xlist,vlist,-h');
            colormap gray(256);
        else
            plot(x, v, 'k.', 'MarkerSize', 12);
        end
        plot(t,-nablaf(t), 'r', 'LineWidth', 2);
        plot([-1 1],[0 0], 'r.', 'MarkerSize', 25);
        axis equal;
        axis([-xmax, xmax, -vmax,vmax]);
        box on;
        set(gca, 'XTick', [], 'YTick', []);
        drawnow;
        mysaveas(k);        
        k = k+1;
        
    end
    v = v + dt/2 * W(x,v);
    x = x + dt*v;
    v = v + dt/2 * W(x,v);
    xvg(:,end+1) = x;
end
