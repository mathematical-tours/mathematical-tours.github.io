rep = 'results/';
if not(exist(rep))
    mkdir(rep);
end

N = 25;

t = linspace(-1,1,N)';

c = [t,t.^2,t.^3];

clf
plot3(c(:,1), c(:,2), c(:,3), 'LineWidth', 2);
axis tight; axis equal;

K = convhull(c(:,1), c(:,2), c(:,3));

clf; hold on;
trisurf(K,c(:,1), c(:,2), c(:,3), t);
view(-30,10);
colormap jet(256);
plot3(c(:,1), c(:,2), c(:,3), 'k', 'LineWidth', 3);
% shading interp;
box on;