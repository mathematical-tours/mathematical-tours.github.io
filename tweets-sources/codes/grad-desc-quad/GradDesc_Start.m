%%
% Gradient descent on a quadratic function, testing anisotropy

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

dotp = @(u,v)sum(u(:).*v(:));

% x^2+a*y^2
a = 2;

f = @(x)(100*(x(2) - x(1)^2)^2 + (1 - x(1))^2);
fxy = @(x,y)100*(y - x.^2).^2 + (1 - x).^2;
nablaf = @(x)[-400*(x(2) - x(1)^2)*x(1) - 2*(1 - x(1)); ...
            200*(x(2) - x(1)^2)];
        
gmode = 'search';
gmode = 'large';
%%
gmode = 'low';
gmode = 'bb';
niter = 200;

q = 50; 
a = .25;

N = 256;
tx = linspace(-2,2,N);
ty = linspace(-1,3,N);
[Y,X] = meshgrid(ty,tx);
    
for it=1:q
    
    t = (it-1)/q*2*pi;
    x = [1.8*cos(t); 1 + 1.8*sin(t)];
    
    F = fxy(X,Y); 
    F = F-min(F(:));
    F = F/max(F(:)); 
    
    
    % initial point
    U = []; 
    g = []; 
    for i=1:niter
        U(:,end+1) = x; 
        % min_t (x-t*gx)^2+a*(y-t*gy)^2
        %    (x-t*x)^2+a*(y-t*a*y)^2
        %    x^2*(1-t)^2+a*y^2*(1-a*t)^2
        %    x^2*(1-t)+y^2*a^2*(1-a*t)=0
        %   t*(x^2+a^3*y^2) = x^2+a^2*y^2
        %   t = (x^2+a^2*y^2)/(x^2+a^3*y^2)
        g1 = g; 
        g = nablaf(x); 
        switch gmode
            case 'search'
                r = (x^2+a^2*y^2)/(x^2+a^3*y^2);
            case 'low'
                r = .7;
                r = .001/4;
            case 'large'
                r = 1.5;
            case 'bb'
                if isempty(g1)
                    r = .5;
                    r = .001/4;
                else
                    r = abs(dotp(x-x1,g-g1))/norm(g-g1)^2;
                end
        end
        x1 = x; 
        x = x - r*nablaf(x);
    end
    
   
    
    s = (it-1)/(q-1);
    s = 0; 
    nc = 15;
    m = linspace(0,1,nc-1)';
    CM = m*[s 0 1-s] + (1-m)*[1 1 1];

    gamma = .35;
    clf; hold on;    
    imagesc(tx,ty,(F').^gamma);
    contour(tx,ty,(F').^gamma,linspace(0,1,nc), 'color', 'k');
    colormap(CM);
    % colormap(parula(256));
   
    plot(U(1,:),U(2,:), 'k.-', 'LineWidth', 2, 'MarkerSize', 23);
    plot(1,1,'r.', 'MarkerSize', 25);
    axis equal; axis([min(tx) max(tx) min(ty) max(ty)]); axis off;
    drawnow;
    
    
    
    saveas(gcf, [rep 'iter-' znum2str(it,2) '.png'], 'png');
    
end
