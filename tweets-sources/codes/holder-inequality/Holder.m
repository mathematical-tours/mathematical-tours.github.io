%%
% Display holder inenquality in 2D.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


n = 256*2;
x = linspace(-1,1,n);
y = linspace(-1,1,n);
[Y,X] = meshgrid(y,x);

p = 1.1;
% 1/p+1/p=1
q = p/(p-1);

x0 = .6; y0 = .9;

N = @(x,y,p)( abs(x).^p+abs(y).^p ).^(1/p);

nanim=60;
for it=1:nanim
    t = (it-1)/(nanim);
    
    x0 = .9*cos(2*pi*t); y0 = .9*sin(2*pi*t);
    
    U = 1 - N(X*x0,Y*y0,1) ./ ( N(X,Y,p) .* N(x0,y0,q) );    
    % U = 1 - (X*x0+Y*y0) ./ ( N(X,Y,p) .* N(x0,y0,q) );
    U = U/max(U(:));
    r = 15; % #levellines   
       
    clf; hold on;
    imagesc(x,y,U');
    if 0
	contour(x,y,U',linspace(0,1,r), 'k');
    colormap(parula(r-1));
    else
    colormap(parula(256));
    end
    caxis([0 1]);
    plot(x0,y0,'.', 'Color', 'r', 'MarkerSize', 20);
    axis equal; axis off;
    drawnow;
    mysaveas(it);
    
end