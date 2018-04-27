function PlotSVM(X,w,b,pi)

col = {[0 0 1] [1 0 0]};
r = .3;
LL = @(c)r*c + (1-r)*[1 1 1]; % lighten a color

ms = 10;

n = [size(X{1},2), size(X{2},2)];
Pi = { pi(1:n(1)), pi(n(1)+1:end) };

% display solution
q = 256; 
u = linspace(0,1,q);
[V,U] = meshgrid(u,u);
R = U*w(1)+V*w(2) + b;
%
rho = 1;
% clf; 
hold on;
imagesc(u,u,-(abs(R')<rho));
contour(u,u,R', [-1 1]*rho, 'k', 'LineWidth', 2);
contour(u,u,R', [0 0], 'k--', 'LineWidth', 2);
colormap(gray(256));
caxis([-3 0]);
axis image; axis equal; axis on; box on;  set(gca, 'XTick', [], 'YTick', []); 
%
for i=1:2
    I = find(abs(Pi{i})>1e-5); %  support vectors
    J = find(abs(Pi{i})<1e-5);
    plot(X{i}(1,J), X{i}(2,J), 'o', 'MarkerFaceColor', LL(col{i}), 'MarkerSize', ms, 'MarkerEdgeColor', col{i});
    plot(X{i}(1,I), X{i}(2,I), 'o', 'MarkerFaceColor', col{i}, 'MarkerSize', ms*1.2, 'MarkerEdgeColor', 'k', 'LineWidth', 2);
end


end