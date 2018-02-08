
rep = '../results/lp-balls/';
[~,~] = mkdir(rep);


q = 1024;
a = 1.05;
t = linspace(-a,a,q);
[Y,X] = meshgrid(t,t);

plist = [0.1:.1:4 50];
plist = [.5 .75 1 1.5 2 4 Inf];
vmax = 1;

for i=1:length(plist)
    p = plist(i);
    if p~=Inf
        F = ( abs(X).^p + abs(Y).^p).^(1/p);
    else
        F = max(abs(X),abs(Y));
    end
    F(F>vmax)=NaN;
    % draw
    r = 15;
    v = linspace(0,vmax,r+2);
    clf; hold on;
    imagesc(t,t,F');
    contour(t,t,F',v, 'k');
    CM = parula(r+1);
    colormap(CM);
    caxis([0,vmax]);
    axis equal; axis off;
    drawnow;
    saveas(gcf, [rep 'ball-' num2str(p) '.png'], 'png');
end