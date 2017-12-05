%%
% Simple example of gradient flow in metric space.


rep = '../results/gradflow-metric/';
[~,~] = mkdir(rep);

f = @(x,y)x.^2 + y.^2;
Gx = @(x,y)2*x;
Gy = @(x,y)2*y;

% background image
N = 256;
tx = linspace(-.05,1,N);
ty = linspace(-.08,.5,N);
[Y,X] = meshgrid(ty,tx);
F = f(X,Y);
clf; hold on;
imagesc(tx,ty,F');
contour(tx,ty,F', 15, 'k');

% for the implicit stepping
N1 = 2048;
tx1 = linspace(min(tx),max(tx),N1);
ty1 = linspace(min(ty),max(ty),N1);
[Y1,X1] = meshgrid(ty1,tx1);


plist = [1 1];
plist = [1 1.3 1.5 1.8 2 3 5 10 Inf];


gmode = 'explicit';
gmode = 'implicit';

tau = .25; niter = 70;
    
lgd = {};
for ip=1:length(plist)
    p=plist(ip);
    
    % conjugate exponent
    q = conjexp(p);
    ri = @(r,p)sign(r).*abs(r).^(1/(p-1));
    
    
    
    grad = @(gx,gy)ri([gx, gy],p);
    if p==1
        grad = @(gx,gy)[gx.*(gx==max(gx,gy)), gy.*(gy==max(gx,gy))];
    elseif p==Inf
        grad = @(gx,gy)[sign(gx),sign(gy)];
    end
    grad = @(gx,gy)grad(gx,gy)*norm([gx,gy],q)/norm(grad(gx,gy),q);
    
    
    % initial point
    x = .95; y = .45;
    U = []; V = [];
    % do the descent
    for i=1:niter
        U(end+1) = x; V(end+1) = y;
        switch gmode
            case 'explicit'
                g = grad(Gx(x,y), Gy(x,y));
                x = x - tau*g(1);
                y = y - tau*g(2);
            case 'implicit'
                E = ( abs(X1-x).^p + abs(Y1-y).^p ).^(2/p) + tau * ( X1.^2 + Y1.^2 );
                if p==Inf
                    E = max(abs(X1-x),abs(Y1-y)).^2 + tau * ( X1.^2 + Y1.^2 );
                end
                [tmp,i] = min(E(:));
                x = X1(i); y = Y1(i);
        end
    end
    
    lgd{end+1} = ['p=' num2str(p)];
    c = (ip-1)/(length(plist)-1);
    plot(U,V, '.-', 'Color', [c 0 1-c], 'LineWidth', 2, 'MarkerSize', 23);
    colormap parula(256);
    axis image; axis off;
    axis([min(tx) max(tx) min(ty) max(ty)]);
    drawnow;
    
end

saveas(gcf, [rep 'gradflow-' gmode '.eps'], 'epsc');