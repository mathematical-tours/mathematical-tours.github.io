addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


n = 256;
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);

q = 75;
plist = linspace(-2,7,q);

    r = 15; % #levellines
    
for it=1:q
    p = plist(it);
    M = ( (X.^p + Y.^p)/2 ).^(1/p);
    %
    s = (it-1)/(q-1);
    m = linspace(0,1,r-1)';
	CM = m*[s 0 1-s] + (1-m)*[1 1 1];
    %
    clf; 
    hold on;
    imagesc(t,t,M);
    contour(t,t,M,linspace(0,1,r), 'k');
    colormap(CM);
    caxis([0 1]);
    axis equal; 
    set(gca, 'XTick', [], 'YTick', []);
    box on;
    drawnow
    mysaveas('anim', it);
end