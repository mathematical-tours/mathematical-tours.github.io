%%
% display various perspective functions
%    psi(x,y) := phi(x/y)*y

name = 'hellinger';
name = 'tv';
name = 'entropy';

addpath('../toolbox/');
rep = MkResRep();



switch name
    case 'entropy'
        phi = @(r)r.*log(r)-r+1;
    case 'tv'
        phi = @(r)abs(r-1);
    case 'hellinger'
        phi = @(r)abs(sqrt(r)-1).^2;
end



r = linspace(0,4,100);
clf;
plot(r, [r.*log(r)-r+1; abs(r-1); abs(sqrt(r)-1).^2], 'LineWidth', 2);
axis tight;
set(gca, 'PlotBoxAspectRatio', [1 1/3 1], 'FontSize', 20);
% set(gca, 'XTick', [0 1], 'YTick', [0 1)]);
saveas(gcf, [rep name '-curves.eps'], 'epsc');
legend('entropy', 'tv', 'hellinger');


psi = @(x,y)phi(x./y).*y;

% f(s*X)=s*f(X)

% y^2/x

x = linspace(1e-4,1,401);
y = linspace(1e-4,1,403);
[Y,X] = meshgrid(y,x);
F = psi(X,Y);

% F(100,:) = 1000;

vmax = min(1,max(F(:)));

% display quantized colormap
r = 15; % #levellines
clf; hold on;
imagesc(x,y,F');
contour(x,y,F',linspace(0,vmax,r), 'k');
plot(x,x, 'k');
colormap(parula(r-1));
caxis([0 vmax]);
axis image; axis off;
saveas(gcf, [rep name '-psi.png'], 'png');