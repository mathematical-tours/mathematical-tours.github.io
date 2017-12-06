% just quadratic functions

rep = 'results/';
[~,~] = mkdir(rep);

% quad matrix
t = .3;
mu = .15;
U = [cos(t) sin(t);-sin(t) cos(t)];
A = U*diag([1,mu])*U';
Ai = inv(A);

% display grid
q = 256;
t = linspace(-1,1,q);
[Y,X] = meshgrid(t,t);


% arrow grid
r = 10;
t1 = linspace(-.9,.9,r);
[Yr,Xr] = meshgrid(t1,t1);

Mlist = {A Ai};
CM = { autumn(256) winter(256) };
for k = 1:length(Mlist)
    B = Mlist{k};
    
    a = 1/2;
    b = -1/4;
    F  = 1/2 * ( B(1,1)*X.^2+B(2,2)*Y.^2+2*B(1,2).*X.*Y ) + a*X + b*Y;
    
    u= B(1,1)*Xr + B(1,2)*Yr - Xr + a;
    v= B(2,1)*Xr + B(2,2)*Yr - Yr + b;
    
    clf; hold on;
    imagesc(t,t,-F');
    contour(t,t,F', 15, 'k');
    colormap(CM{k});
    caxis([min(-F(:)) max(-F(:))]);
    axis off;
    
    rho = 1;
    quiver(Xr,Yr,rho*u,rho*v, rho, 'k', 'MaxHeadSize', .05);
    
end;