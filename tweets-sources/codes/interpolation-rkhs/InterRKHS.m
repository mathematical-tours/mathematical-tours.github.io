%%
% Interpolation using RKHS kernels.

addpath('../toolbox/');
rep = MkResRep();
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


kernel = 'cubic';
kernel = 'gaussian';
kernel = 'mutiquadric';
kernel = 'invquadric';
kernel = 'matern';

switch kernel
    case 'gaussian'
        phi = @(s)exp(-s.^2/2);
    case 'laplacian'
        phi = @(s)exp(-abs(s));
    case 'matern'
        phi = @(s)exp(-s).*(1+s);
    case 'mutiquadric'
        phi = @(s)sqrt(1+s.^2);
    case 'invquadric'
        phi = @(s)1./(1+s.^2);
    case 'polyharmonic'
        phi = @(s)s.^3;
end
K = @(X,Y)phi(distmat(X,Y));

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

% D0 = distmat(x', x');
% D = distmat(X(:)', x');

slist = [.02 .03 .05 .075 .1 .2];

clf; hold on;
for i=1:length(slist);
    s = slist(i);
    t=(i-1)/(length(slist)-1);
    Q = K(x'/s,x'/s) \ f; % phi(D0/s) \ f; 
    F = K(X(:)'/s,x'/s) * Q; % phi(D/s)*Q;
    %
    plot(X,F, 'LineWidth', 2, 'Color', [t 0 1-t]);
end
plot(x,f, '.k', 'MarkerSize', 25);
axis tight;
box on;
SetAR(2/3);
axis([0,1,-.05,1.05]);
set(gca, 'XTick', [], 'YTick', []);
saveas(gcf, [rep kernel '-1d.eps'], 'epsc');


%%
% In 2D.


n = 256; % grid

if not(exist('Z'))
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
Z = [real(Z); imag(Z)];
end
m = size(Z,2);

% evaluation point
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
D0 = distmat(Z, Z);
D = distmat([X(:), Y(:)]', Z);

% discrete values
rand('state', 1346);
f = rescale(rand(m,1), -1, 1); 

for i=1:length(slist)
    s = slist(i);
    % solve
    Q = phi(D0/s) \ f; 
    F = phi(D/s)*Q;
    F = reshape(F, [n n]);
    % plot
    r = 30; % #levellines
    clf; hold on;
    imagesc(t,t,F');
    contour(t,t,F',linspace(min(F(:)),max(F(:)),r), 'k');
    colormap(parula(r-1));
    caxis([min(F(:)),max(F(:))]);
    plot(Z(1,:), Z(2,:), 'r.', 'MarkerSize', 20);
    axis image; axis off;
    saveas(gcf, [rep kernel '-2d-q' num2str(s) '.png'], 'png'); 
end