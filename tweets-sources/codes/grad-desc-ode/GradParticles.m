
addpath('../toolbox/');
rep = MkResRep();

rho = 0.15;
peak =  @(x,y)3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ...
   - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ...
   - 1/3*exp(-(x+1).^2 - y.^2) + rho*x.^2 + rho*y.^2;

h = 1e-6;
peakg = @(x,y)cat(3, (peak(x+h,y)-peak(x,y))/h, (peak(x,y+h)-peak(x,y))/h);

a = 2.5;
n = 501;
t = linspace(-a,a,n);
[Y,X] = meshgrid(t,t);
R = peak(X,Y);

% points for gradient display
q = 20;
t1 = linspace(-a,a,q+2); t1 = t1(2:end-1);
[Y1,X1] = meshgrid(t1,t1);
G = peakg(X1,Y1);

u = min(R(:)); v = max(R(:));
r = 25; % #levellines
clf; hold on;
imagesc(t,t,R');
contour(t,t,R',linspace(u,v,r), 'k');
colormap(parula(r-1));
caxis([u v]);
axis image; axis off;
% display arrows
quiver(X1,Y1,G(:,:,1), G(:,:,2), 1.5, 'k');
saveas(gcf, [rep 'gradient.png']);

% integrate ODE in time
b = a*.9;
%
m = 200; % #particules
P0 = rand(m,1,2)*2*b-b;
%
q = 15;
q = 20;
m = q*q;
t1 = linspace(-b,b,q); 
[Y1,X1] = meshgrid(t1,t1);
P0 = reshape( [X1(:),Y1(:)], [m 1 2] );

P0 = randn(30000,1,2);



tau = .01;
niter = 200;

tau = .025;
niter = 80;
%
P = P0;
for i=1:niter
    g = peakg(P(:,end,1), P(:,end,2));
    % normalize
    % g = g ./ repmat(sqrt(sum(g.^2,3)), [1 1 2]);
    g = g ./ max(max(sqrt(sum(g.^2,3))));
    %
    P(:,end+1,:) = P(:,end,:) - tau*g;    
end

ndisp = min(3000,size(P0,1));
for k=1:niter
    c = (k-1)/(niter-1);
    % display
    clf; hold on;
    imagesc(t,t,R');
    contour(t,t,R',linspace(u,v,r), 'k');
    colormap(parula(r-1));
    caxis([u v]);
    axis image; axis off;
    % display trajectories
    lw = 1;
%    for i=1:m
%        PlotTrajectory( P(i,1:k,1), P(i,1:k,2), lw, [0 0 1], [c 0 1-c] );
%    end
    plot( P(1:ndisp,k,1), P(1:ndisp,k,2), '.', 'color', [c 0 1-c], 'MarkerSize', 15);
    axis equal;
    axis([-a a -a a]); 
    saveas(gcf, [rep 'points-' znum2str(k,2) '.png']);
    %    
    % render density
    sc = 10;
    h = parzen2d( .5*P(:,k,1)/b+.5, .5*P(:,k,2)/b+.5,256,sc); 
    h = h/max(h(:));
    m = linspace(0,1,r-1)';
    CM = m*[c 0 1-c] + (1-m)*[1 1 1];
    clf; hold on;
    tt = linspace(0,1,256);
    imagesc(tt,tt,h');
    contour(tt,tt,h',linspace(0,1,r), 'Color',[c 0 1-c]);
    colormap(CM);
    caxis([0 1]); box on;
    axis equal; set(gca, 'XTick', [], 'YTick', []);
    % drawnow; 
    % saveas(gcf, [rep 'density-' znum2str(k,2) '.png']);
end