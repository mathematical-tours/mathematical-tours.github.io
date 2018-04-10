%%
% Interpolation using natural neighbors.


addpath('../toolbox/');
rep = MkResRep();
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

m = 12;
n = 300; % grid


rand('state', 1234);
Z = rand(2,m);

% click selection
it = 0;
clf; hold on;
Z = [];
while true
    axis([0 1 0 1]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    Z(end+1) = a+1i*b;
end
m = length(Z);
Z = [real(Z); imag(Z)];



% evaluation point
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
Dz0 = distmat([X(:), Y(:)]', Z);
[Dz,Iz] = min(Dz0,[], 2);

% compute the weights
W = zeros(n*n,1);
for i=1:n*n
    progressbar(i,n*n);
    x = [X(i); Y(i)];
    %
    Dx = distmat([X(:), Y(:)]', x);
    D = min(Dx,Dz);
    J = find(D==Dx);  % new voronoi cells
    for k=1:m
        W(i,k) = sum(Iz(J)==k)/length(J);
    end
end

% use them to interpolate values


f = rand(m,1);
F = reshape(W*f, [n n]);

% plot
r = 25; % #levellines
clf; hold on;
imagesc(t,t,F');
contour(t,t,F',linspace(min(F(:)),max(F(:)),r), 'k');
colormap(parula(r-1));
caxis([min(F(:)),max(F(:))]);
plot(Z(1,:), Z(2,:), 'r.', 'MarkerSize', 25);
axis image; axis off;
saveas(gcf, [rep 'natural-m' num2str(m) '.png'], 'png');

% display the weights
for k=1:m
    w = reshape(W(:,k), [n n]);
    clf; hold on;
    imagesc(t,t,w');
    contour(t,t,w',linspace(min(w(:)),max(w(:)),r), 'k');
    colormap(parula(r-1));
    caxis([min(w(:)),max(w(:))]);
    plot(Z(1,:), Z(2,:), 'r.', 'MarkerSize', 20);
    plot(Z(1,k), Z(2,k), 'b.', 'MarkerSize', 30);
    axis image; axis off;
    saveas(gcf, [rep 'natural-m' num2str(m) '-' num2str(k) '.png'], 'png');
end
