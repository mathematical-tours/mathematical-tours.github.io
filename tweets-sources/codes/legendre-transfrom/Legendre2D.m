%%
% Simple legendre transforms.


addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name, it)saveas(gcf, [rep  'a-' name '-' znum2str(it,3) '.png']);


nb = 70;

P = 1 + 1e-3 + 10 * linspace(0,1, nb).^2;
Q = P./(P-1); % 1/q = 1-1/p ==> q = P./(P-1);

n = 256;
x = linspace(-3,3,n);
[Y,X] = meshgrid(x,x);
R = sqrt(X.^2+Y.^2);


% 

for it=1:nb
    p = P(it);
    t = (it-1)/(nb-1);
    F = R.^p./p;
    clf; 
    surf(x,x, F, rescale(min(F,3)) ); % , 'LineWidth', 3, 'Color', [t 0 1-t] );
    view(3); shading interp;
    axis([min(x) max(x) min(x) max(x) 0 3]); box on; 
    set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    camlight;
    drawnow;
    mysaveas('primal', it);
end

for it=1:nb
    q = Q(it);
    t = (it-1)/(nb-1);
    F = R.^q./q;
    clf; 
    surf(x,x, F, rescale(min(F,3)) ); % , 'LineWidth', 3, 'Color', [t 0 1-t] );
    view(3); shading interp;
    axis([min(x) max(x) min(x) max(x) 0 3]); box on; 
    set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    camlight;
    drawnow;
    mysaveas('dual', it);
end