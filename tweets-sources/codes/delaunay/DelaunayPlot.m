%%
% Plot of Delaunay triangulation

rep = '../results/delaunay/';
[~,~] = mkdir(rep);
addpath('../toolbox/');

if not(exist('test'))
    test = 0;
end
test = test+1;


% input points
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
% add corner
Z(:,end+1:end+4) = [[1 1 0 0];[1 0 0 1]];
Zc = Z(1,:) + 1i*Z(2,:);

% compute delaunay
T = delaunay(Z(1,:),Z(2,:));
[VX,VY] = voronoi(Z(1,:),Z(2,:));
V = VX + 1i*VY;
% circum centers
DT = delaunayTriangulation(Z');
W = circumcenter(DT);
Wc = W(:,1)' + 1i*W(:,2)';
% circum radius
r = abs(Wc - Zc(T(:,1)));

t = linspace(0,1,100);
c = exp(2i*pi*t);

%%
% Display Delaunay

clf; hold on;
triplot( T, Z(1,:), Z(2,:), 'r.-', 'LineWidth', 2, 'MarkerSize', 25 );
plot(V, '--', 'LineWidth', 1, 'color', 'b');
axis equal; axis([0 1 0 1]); 
axis off;
% circum circles
if size(Z,2)<=10
plot(Wc, 'k.', 'MarkerSize', 20);
for i=1:length(W)
    plot( Wc(i)+r(i)*c, 'k', 'LineWidth', 2 );
end
end
saveas(gcf, [rep num2str(test) '-delaunay.eps'], 'epsc');


%%
% Display Voronoi

clf; hold on;
triplot( T, Z(1,:), Z(2,:), 'r.--', 'LineWidth', 1, 'MarkerSize', 25 );
plot(V, 'b-', 'LineWidth', 2);
axis equal; axis([0 1 0 1]); 
axis off;
saveas(gcf, [rep num2str(test) '-voronoi.eps'], 'epsc');