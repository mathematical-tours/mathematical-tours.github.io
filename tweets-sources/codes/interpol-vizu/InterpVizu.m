%%
% Display in 3D interpolation

addpath('../toolbox/');
rep = MkResRep();

n = 5; 
k = 6;

x  = 1/2:n-1/2; %   linspace(0,1,n);
xi = linspace(0,n,k*n);

x  = linspace(0,1,n);
xi = linspace(0,1,k*n);

rand('state', 12);
f0 = rescale( rand(n) + x );
f1 = rescale( rand(n) + x' );

meth = 'cubic';
meth = 'linear';
meth = 'nearest';

q = 50;
for it=1:q
    s = (it-1)/(q-1);
    switch meth
            spapi( optknt(x,k), x, y );
        otherwise            
            fi = interp2(x,x',(1-s)*f0+s*f1,xi,xi',meth);    
    end
    clf;
    surf(xi,xi,fi);
    axis([0 1 0 1 0 1]);
    colormap parula(256);
    caxis([0 1]);
    axis off;
    drawnow;
    % saveas(gcf, [rep meth '-' znum2str(it,2) '.png'] );
end