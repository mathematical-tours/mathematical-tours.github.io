%%
% Interpolation using Shepard weights.


addpath('../toolbox/');
rep = MkResRep();
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


%% 
% In 1D.

% points and values
m = 6;  % #points
n = 256*2; % grid

rand('state', 123);
x = rand(m,1)+.5;
x = rescale(cumsum(x), .1, .9); 
f = rescale(rand(m,1));


X = linspace(0,1,n);
D = distmat(X(:)', x');

qlist = [1 1.5 2 3 6 100];
clf; hold on;
for i=1:length(qlist);
    q = qlist(i);
    t=(i-1)/(length(qlist)-1);
    %
    phi = @(r)1./(r+1e-10).^q;
    F = ( phi(D)*f ) ./ sum( phi(D), 2 );
    %
    plot(X,F, 'LineWidth', 2, 'Color', [t 0 1-t]);
end
plot(x,f, '.k', 'MarkerSize', 25);
axis tight;
box on;
SetAR(2/3);
axis([0,1,-.05,1.05]);
set(gca, 'XTick', [], 'YTick', []);
saveas(gcf, [rep 'shepard-1d.eps'], 'epsc');


%%
% In 2D.

m = 200;  % #points
m = 12;
n = 256; % grid

r = 3; % BB size
rand('state', 1234);
xy = (2*rand(2,m)-1)*r;

% evaluation point
t = linspace(-r,r,n);
[Y,X] = meshgrid(t,t);
D = distmat([X(:), Y(:)]', xy);

% discrete values
f = peaks(xy(1,:), xy(2,:)); f = f(:);
f = rand(m,1); 

qlist = [1 1.5 2 3 6 100];
for i=1:length(qlist)
    q = qlist(i);
    % solve
    phi = @(r)1./(r+1e-10).^q;
    F = ( phi(D)*f ) ./ sum( phi(D), 2 );
    F = reshape(F, [n n]);
    % plot
    r = 30; % #levellines
    clf; hold on;
    imagesc(t,t,F');
    contour(t,t,F',linspace(min(F(:)),max(F(:)),r), 'k');
    colormap(parula(r-1));
    caxis([min(F(:)),max(F(:))]);
    plot(xy(1,:), xy(2,:), 'r.', 'MarkerSize', 20);
    axis image; axis off;
    saveas(gcf, [rep 'shepard-2d-q' num2str(q) '.png'], 'png'); 
end
