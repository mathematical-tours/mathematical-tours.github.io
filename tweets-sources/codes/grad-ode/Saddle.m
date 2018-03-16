%%
% Display gradient flow on saddle points

rep = '../results/saddles/';
[~,~] = mkdir(rep);

a = 3; b = 3;
a = 2; b = 3;
a = 2; b = 2;
a = 4; b = 3;
a = 1; b = 3;

f = @(x,y)x.^a/a - y.^b/b;
h = 1e-6;
fG = @(x,y)cat(3, x.^(a-1), -y.^(b-1));

L = 1;
n = 501;
t = linspace(-L,L,n);
[Y,X] = meshgrid(t,t);
R = f(X,Y);

% points for gradient display
q = 20;
t1 = linspace(-L,L,q+2); t1 = t1(2:end-1);
[Y1,X1] = meshgrid(t1,t1);
G = fG(X1,Y1);

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
La = L*.9;
%
q = 10;
m = q*q;
t1 = linspace(-La,La,q); 
[Y1,X1] = meshgrid(t1,t1);
P0 = reshape( [X1(:),Y1(:)], [m 1 2] );

m = 52; 
t1 = (0:m-1)'/m*2*pi;
P0= reshape( .75*[cos(t1),sin(t1)], [m 1 2] );


% gradient descent
tau = .2;
P = P0;
niter = 1000;
for i=1:niter
    g = fG(P(:,end,1), P(:,end,2));
    P(:,end+1,:) = P(:,end,:) - tau*g;    
end

clf; hold on;
imagesc(t,t,R');
contour(t,t,R',linspace(u,v,r), 'k');
colormap(parula(r-1));
caxis([u v]);
axis image; axis off;

% display trajectories
for i=1:m
    A = squeeze(P(i,:,:));
    I = find(max(abs(A),[],2)<=L*1.2);
    PlotTrajectory( P(i,I,1), P(i,I,2) );
%    plot(P(i,I,1), P(i,I,2), 'k-', 'LineWidth', 2);
    drawnow;
end
axis([-L L -L L]); axis equal;




saveas(gcf, [rep 'flow-' num2str(a) '-' num2str(b) '.png']);

