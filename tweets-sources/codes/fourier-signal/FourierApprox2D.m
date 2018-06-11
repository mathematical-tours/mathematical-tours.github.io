%%
% Approximation of signals using fft.

addpath('../toolbox/');
rep = MkResRep('image');

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

n = 512;
name = 'hibiscus';
f = rescale(sum(load_image(name, n),3));
[fS,I] = sort(f(:)); f(I) = linspace(0,1,length(f(:)));

m = [0:n/2,-n/2+1:-1]';
[Y,X] = meshgrid(m,m);
R = sqrt(X.^2+Y.^2);

q = 70;
plist = 1 + (n/12-1)*linspace(0,1,q);

u = linspace(0,1,n);

for i=1:q
    t = (i-1)/(q-1);
    p = plist(i);
    M = R<=p;
    f1 = real(ifft2( fft2(f).*M ) );
    % display
    r = 12; % #levellines
    clf; hold on;
    imagesc(u,u,f1');
    contour(u,u,f1',linspace(0,1,r), 'k');
    colormap(gray(r-1));
    caxis([0 1]);
    plot([0 1 1 0 0 1], [0 0 1 1 0 0], 'color', [t 0 1-t], 'LineWidth', 5)
    axis image; axis off;
    drawnow;
    saveas(gcf, [rep name '-' znum2str(i,2) '.png']);
end

% AutoCrop(rep, name);