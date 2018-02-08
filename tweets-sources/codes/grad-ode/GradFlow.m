%%
% Display gradient flow

rep = '../results/grad-flow/';
[~,~] = mkdir(rep);

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
m = q*q;
t1 = linspace(-b,b,q); 
[Y1,X1] = meshgrid(t1,t1);
P0 = reshape( [X1(:),Y1(:)], [m 1 2] );



tau = .01;
P = P0;
niter = 200;
for i=1:niter
    g = peakg(P(:,end,1), P(:,end,2));
    P(:,end+1,:) = P(:,end,:) - tau*g;    
end

clf; hold on;
imagesc(t,t,R');
contour(t,t,R',linspace(u,v,r), 'k');
colormap(parula(r-1));
caxis([u v]);
axis image; axis off;
% display trajectories
% plot(P(:,:,1)', P(:,:,2)', 'k');
for i=1:m
    PlotTrajectory( P(i,:,1), P(i,:,2) );
end
axis([-a a -a a]); axis equal;
saveas(gcf, [rep 'flow.png']);

