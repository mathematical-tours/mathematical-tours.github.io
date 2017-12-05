%%
% Mirror descent

rep = 'results/';
[~,~] = mkdir(rep);


N = 256;
tx = linspace(1e-3,1,N);
ty = linspace(1e-3,1,N);
[Y,X] = meshgrid(ty,tx);

breglist = {'kl' 'sqrt' 'log' 'eucl'};


x = .6; y = .6;
    
for it = 1:length(breglist)
    
    bregmode = breglist{it};
    
    [R,Ri,phi,D] = load_bregman(bregmode);
    
    F = D(X,x*ones(size(X))) + D(Y,y*ones(size(X)));
    F = sqrt(F);
    clf; hold on;
    imagesc(tx,ty,-F');
    contour(tx,ty,F', 15, 'k');
    colormap parula(256);
    caxis([min(-F(:)),max(-F(:))]);
    axis off; axis equal;
    plot(x,y,'r.', 'MarkerSize', 30);
    saveas(gcf, [rep 'bregman-' bregmode '.png'], 'png');
    
end





