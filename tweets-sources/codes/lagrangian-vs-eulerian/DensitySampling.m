%%
% Contrast Eulerian vs Lagrangian discretization.

addpath('../toolbox/');
rep = MkResRep();

% load initial density
n = 512;

t = (0:n-1)'/n;
[Y,X] = meshgrid(t,t);

f = peaks(n);
f = abs(f);
f = f/sum(f(:));

%%
% display quantized colormap

r = 15; % #levellines
clf; hold on;
imagesc(t,t,1-f'/max(f(:)));
contour(t,t,f'/max(f(:)),linspace(0,1,r), 'k');
colormap(gray(r-1));
caxis([0 1]);
axis image; 
set(gca, 'XTick', [], 'YTick', [], 'box', 'on');
saveas(gcf, [rep 'density.png']);

%% 
% Eulerian.

q = 60; % 
Nlist = 1:q;
for i=1:q
    g = rescale( ApproxInterp(f,Nlist(i)) );
    c = (i-1)/(q-1); col = [c 0 1-c];
    G = zeros(n,n,3);
    for k=1:3
        G(:,:,k) = g*col(k) + (1-g)*1;
    end
    clf; imagesc(G); axis off; axis image; drawnow;
    imwrite(G, [rep 'eulerian-' znum2str(i,2), '.png']);
end

%% 
% Lagrangian.

q = 60; % 
Nlist = round(linspace(4,3000,q));
I = randsample(n*n,max(Nlist),true,f(:));

t = (0:n-1)'/n;
[Y,X] = meshgrid(t,t);
ms = 20;
for i=1:q
    c = (i-1)/(q-1); col = [c 0 1-c];
    m = Nlist(i);
    x = X(I(1:m)); y = Y(I(1:m));
    clf; hold on;
    a = ones(m,1)*35; alp = .2;
    scatter(x,y,a,'MarkerFaceColor',col, ... %'MarkerEdgeColor',col,...
        'MarkerFaceAlpha',alp,'MarkerEdgeAlpha',0);    
    axis equal;
    axis([0 1 0 1]);  box on;
    set(gca, 'XTick', [], 'YTick', [], 'box', 'on');
    axis on; axis ij,
    drawnow;
    saveas(gcf, [rep 'lagrangian-' znum2str(i,2) '.png']);
end
% AutoCrop(rep, 'lagrangian-');
