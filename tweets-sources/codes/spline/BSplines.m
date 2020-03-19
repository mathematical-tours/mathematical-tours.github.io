%% 
% Periodic b-splines 

addpath('../toolbox/');
rep = MkResRep();

clf; hold on;
Z = [];
while true
    axis equal; axis([0 1 0 1]);  % axis off;
    [a,b,button] = ginput(1);  
    if button==3
        break;
    end
    plot(a,b, '.', 'MarkerSize', 25); 
    Z(end+1) = a + 1i*b;
    plot(Z, 'r', 'LineWidth', 2);
end
n = size(Z,2);
Z = Z(:);
W = [];
for i=1:n
    axis equal; axis([0 1 0 1]);  % axis off;
    [a,b,button] = ginput(1);  
    plot(a,b, '.', 'MarkerSize', 25); 
    W(end+1) = a + 1i*b;
    plot(W, 'b', 'LineWidth', 2);
end
W = W(:);

p = length(Z);
n = 3; % order of the splines
m = 1024; % for vizualization

t = (0:p-1)/p; % knot
x = (0:m-1)'/m; %linspace(0,1,m)';
B = SplineBasis(x,t,n,p);


% the basis
clf;
plot(x, B, 'LineWidth', 2);
axis([0 1 0 max(B(:))*1.05]);
box on;set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1]);
saveas(gcf, [rep 'basis.png']);
drawnow;

q = 70;
for it=1:q
    s = (it-1)/(q-1);
    %
    X = (1-s)*Z + s*W;
    y = B*X;
    %
    clf; hold on;
    plot(X([1:end 1]), 'k.-', 'LineWidth', 1, 'MarkerSize', 20);
    plot(y, 'r', 'Color', [s 0 1-s], 'LineWidth', 2);
    axis equal; axis([0 1 0 1]);
    box on; set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1]);
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
    drawnow;
end


