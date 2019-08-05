%%
% Display gradient flow


addpath('../toolbox/');
rep = MkResRep('3d');

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


r = 12;
ls = linspace(min(R(:)), max(R(:)), r);
% extract levelsets
[C,h] = contour(t,t,R,ls);
CL = {};
while not(isempty(C))
    m = C(2,1); C(:,1) = [];
    c = C(:,1:m); C(:,1:m) = [];
    c(3,:) = peak(c(2,:),c(1,:));
    CL{end+1} = c;
end

% display
clf; hold on;
surf(t,t,R');
for i=1:length(CL)
    cl = CL{i};
    plot3(cl(2,:), cl(1,:), cl(3,:), 'k', 'LineWidth', 2);
end
shading interp;
view(-65,40);
colormap(parula(r-1));
caxis([min(R(:)), max(R(:))]);
camlight; axis off;


% % integrate ODE in time
% b = a*.9;
% %
% q = 9;
% m = q*q;
% t1 = linspace(-b,b,q); 
% [Y1,X1] = meshgrid(t1,t1);
% P0 = reshape( [X1(:),Y1(:)], [m 1 2] );


% click selection
if not(exist('P0'))
it = 0; P0 = [];
clf; hold on;
imagesc(t,t,R');
while true
    [u,v,button] = ginput(1);
    plot(u,v, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    P0(end+1,1,1) = u;
    P0(end,1,2) = v;
end
end
m = size(P0,1);


tau = .01;
niter = 200;

tau = .025;
niter = 70;
%
P = P0;
for i=1:niter
    g = peakg(P(:,end,1), P(:,end,2));
    % normalize
    g = g ./ repmat(sqrt(sum(g.^2,3)), [1 1 2]);
    %
    P(:,end+1,:) = P(:,end,:) - tau*g;    
end

% delta for each frame in term of rotation
delta_rot = 2;

for k=1:niter
    c = (k-1)/niter;
    
    % display
    clf; hold on;
    surf(t,t,R');
    for i=1:length(CL)
        cl = CL{i};
        plot3(cl(2,:), cl(1,:), cl(3,:), 'k', 'LineWidth', 2);
    end
    shading interp;
    view(-65+delta_rot*(k-1),40);
    colormap(parula(r-1));
    caxis([min(R(:)), max(R(:))]);
    camlight; axis off;
    % plot3( P(:,1,1), P(:,1,2), peak(P(:,1,1), P(:,1,2)), '.', 'color', 'r', 'MarkerSize', 30);
    
    % display trajectories
    lw = 2;
    for i=1:m
        PlotTrajectory3D( P(i,1:k,1), P(i,1:k,2), peak(P(i,1:k,1), P(i,1:k,2)), lw, [0 0 1], [c 0 1-c] );
    end
    axis([-a a -a a min(R(:)) max(R(:))]);
    % plot3( P(:,k,1), P(:,k,2), peak(P(:,k,1), P(:,k,2)), '.', 'color', [c 0 1-c], 'MarkerSize', 25);
    saveas(gcf, [rep 'flow-' znum2str(k,2) '.png']);
end

% AutoCrop(rep, ['flow-']); 
