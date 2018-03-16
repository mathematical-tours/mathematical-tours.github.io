%%
% DC optimization

rep = '../results/dc-optim/';
[~,~] = mkdir(rep);
addpath('../toolbox/');

x = linspace(-2,2,1023);

plot(x,  x.^4/4 - x.^2/2, 'LineWidth', 2 );

u = 1.3;
t = linspace(-u,u,512);
[Y,X] = meshgrid(t,t);


D = @(k)(X.^k + Y.^k)/k;

F = D(4)-D(2);
% FS = D(2)-D(4);


r = 15; % #levellines
clf; hold on;
imagesc(t,t,F);
contour(t,t,F,linspace(min(F(:)),max(F(:)),r), 'k');
colormap(parula(r-1));
caxis([min(F(:)) max(F(:))]);
axis image; axis off;


K = 50; 
C = distinguishable_colors(K);
for k=1:K
    % optimization algorithm
    niter = 20;
    x = [];
    nablaF  = @(x)x;
    cr = @(s)sign(s).*abs(s).^(1/3);
    nablaGs = @(x)cr(real(x)) + 1i*cr(imag(x));
    x = (randn + 1i*randn)/5;
    x = .1 + 1.2i;
    [a,b,button] = ginput(1);
    if button==3
        break;
    end
    x = a+1i*b;
    y = [];
    for i=1:niter
        y(end+1) = nablaF(x(end));
        x(end+1) = nablaGs(y(end));
    end
    plot(x, 'r.-', 'LineWidth', 2, 'MarkerSize', 25, 'Color', C(k,:));
end
saveas(gcf, [rep 'dcoptim.png']);

return;

clf; hold on;
imagesc(t,t,FS);
contour(t,t,FS,linspace(min(FS(:)),max(FS(:)),r), 'k');
colormap(parula(r-1));
caxis([min(FS(:)) max(FS(:))]);
axis image; axis off;
plot(y, 'r.-', 'LineWidth', 2, 'MarkerSize', 20);

