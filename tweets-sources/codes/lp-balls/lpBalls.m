

addpath('../toolbox/');
rep = MkResRep('2d');

n = 1024;
a = 1.05;
t = linspace(-a,a,n);
[Y,X] = meshgrid(t,t);

q = 50; 
plist = linspace(.5,5,q);
vmax = 1;

%% In 2D

for i=1:q
    s = (i-1)/(q-1);
    %
    p = plist(i);
    if p~=Inf
        F = ( abs(X).^p + abs(Y).^p).^(1/p);
    else
        F = max(abs(X),abs(Y));
    end
    %
    imAlpha = ones(n);
    imAlpha(F>vmax) = 0;
    % generate colormap
    r = 15;
    CM = parula(r+1);
    m = linspace(0,1,r+1)';
    CM = m*[s 0 1-s] + (1-m)*[1 1 1];
    % draw
    v = linspace(0,vmax,r+2);
    clf; hold on;
    imagesc(t,t,F', 'AlphaData', imAlpha);
    contour(t,t,F',v, 'k');
    % contour(t,t,F',[vmax vmax], 'LineWidth', 3, 'Color', [s 0 1-s]);    
    colormap(CM);
    caxis([0,vmax]);
    axis equal; axis off;
    drawnow;
    saveas(gcf, [rep '2d-balls-' znum2str(i,2) '.png'], 'png');
end


% AutoCrop(rep, ['2d-']); 


