%%
% Gradient descent on a quadratic function.


addpath('../toolbox/');
rep = MkResRep();

% x^2+a*y^2
a = 2;

gmode = 'large';
gmode = 'low';
gmode = 'search';
niter = 10;

q = 50; 
alist = linspace(.01,5,q);

alist = linspace(1,50,q/2);
alist = .5 + 60*linspace(0,1,q).^2;

N = 256;
tx = linspace(-.3,1,N);
ty = linspace(-.6,.6,N);
[Y,X] = meshgrid(ty,tx);
    
for it=1:q
    
    a = alist(it);
    
    F = ( X.^2+a*Y.^2 )/2;
    
    
    % initial point
    x = .9; y = .3;
    x = .9; y = .3;
    U = []; V = [];
    for i=1:niter
        U(end+1) = x; V(end+1) = y;
        % min_t (x-t*gx)^2+a*(y-t*gy)^2
        %    (x-t*x)^2+a*(y-t*a*y)^2
        %    x^2*(1-t)^2+a*y^2*(1-a*t)^2
        %    x^2*(1-t)+y^2*a^2*(1-a*t)=0
        %   t*(x^2+a^3*y^2) = x^2+a^2*y^2
        %   t = (x^2+a^2*y^2)/(x^2+a^3*y^2)
        switch gmode
            case 'search'
                r = (x^2+a^2*y^2)/(x^2+a^3*y^2);
            case 'low'
                r = .1;
            case 'large'
                r = .52;
        end
        x = x - r*x;
        y = y - r*a*y;
    end
    
    
    plot(log(abs(U+1i*V)));
    
    s = (it-1)/(q-1);
    nc = 15;
    m = linspace(0,1,nc-1)';
    CM = m*[s 0 1-s] + (1-m)*[1 1 1];

    
    clf; hold on;
    imagesc(tx,ty,F'/max(F(:)));
    contour(tx,ty,F'/max(F(:)),linspace(0,1,nc), 'color', [s 0 1-s]);
    plot(U,V, 'k.-', 'LineWidth', 2, 'MarkerSize', 23);
    colormap(CM);
    axis equal; axis([min(tx) max(tx) min(ty) max(ty)]); axis off;
    drawnow;
    
    
    
    saveas(gcf, [rep 'iter-' znum2str(it,2) '.png'], 'png');
    
end
