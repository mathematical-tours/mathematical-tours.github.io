function mydisp(M)

n = size(M,1);
t = linspace(0,1,n);

% display quantized colormap
r = 15; % #levellines
hold on;
imagesc(t,t,M);
contour(t,t,M,linspace(0,1,r), 'k');
colormap(parula(r-1));
caxis([0 1]);
axis image; axis off;


end
