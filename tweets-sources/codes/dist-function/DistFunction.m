%%
% test for computation of distance functions.

rep = '../results/dist-function/';
[~,~] = mkdir(rep);
addpath('../toolbox/');


name = 'curvepoints';
name = 'points';
name = 'curve';

clf; hold on;
Z = [];
while true
    axis([0 1 0 1]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    Z(:,end+1) = [a;b];
end

if strcmp(name, 'curve') || strcmp(name, 'curvepoints')
    m = 0;
    if strcmp(name, 'curvepoints')
        m = 4;
    end
    % perform curve subdivision
    subdivide = @(f,h)cconvol( upsampling(f), h);
    h = [1 4 6 4 1]; % cubic B-spline
    h = 2*h/sum(h(:));
    Z0 = Z(:,end-m+1:end);
    z = Z(1,1:end-m)'+1i*Z(2,1:end-m)';
    for k=1:6
        z = subdivide(z,h);
    end
    Z = [real(z(:))'; imag(z(:))'];
else
    Z0 = Z; Z = [];
end

% grid
n = 256*2;
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
% distance
D = distmat( [Z Z0],[X(:)';Y(:)'] );
D = reshape(min(D,[],1),[n n]);

% draw
r = 15;
clf; hold on;
imagesc(t,t,D');
contour(t,t,D',r, 'k');
colormap(parula(r+1));
%
plot(Z0(1,:), Z0(2,:), 'r.', 'MarkerSize', 25);
if not(isempty(Z))
    plot(Z(1,[1:end 1]), Z(2,[1:end 1]), 'r', 'LineWidth', 2);
end
axis tight; axis equal;
axis off;
drawnow;
saveas(gcf, [rep name '-dist.png'], 'png');

% compute delaunay
Z1 = [Z Z0];
z = Z1(1,:) + 1i*Z1(2,:);
T = delaunay(Z1(1,:),Z1(2,:))';
[VX,VY] = voronoi(Z1(1,:),Z1(2,:));
V = VX + 1i*VY;

% draw
r = 15;
clf; hold on;
imagesc(t,t,D');
% contour(t,t,D',r, 'k--');
colormap(parula(r+1));
%
plot(Z0(1,:), Z0(2,:), 'r.', 'MarkerSize', 25);
if not(isempty(Z))
    plot(Z(1,[1:end 1]), Z(2,[1:end 1]), 'r', 'LineWidth', 2);
end
% medial axis
m = abs(V(1,:)-V(2,:));
if strcmp(name, 'points')
    s = Inf
else
    s = .1; % threshold
end
I = find(m<s);
plot(V(:,I), '-', 'LineWidth', 3, 'color', 'b');
%
axis tight; 
axis([0 1 0 1]); axis equal;
axis([0 1 0 1]);
axis off;
drawnow;
saveas(gcf, [rep name '-skel.png'], 'png');



