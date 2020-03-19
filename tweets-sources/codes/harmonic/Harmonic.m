% harmonic functions 2D

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name, it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


n = 160;

% points indicating constraints
k = 4;
x = rand(k,1)+1i*rand(k,1);
a = (-1).^(1:k);
% speed 
eta = .03;
v = randn(k,1) + 1i*randn(k,1);
v = eta*v./abs(v);

proj = @(u)clamp(real(u),1,n) + 1i * clamp(imag(u),1,n);
swap = @(u)imag(u) + 1i*real(u);

q = 100;
niter = 30000;
f = zeros(n);
f = HarmDiffus(f, proj(ceil(x*n)), a, niter, 300);
niter = 5000;

for it=1:q
    
    X = proj(ceil(x*n));
    f = HarmDiffus(f, X, a, niter, Inf);
    
    
    
    % display quantized colormap
    R = (f+1)/2;
    t = linspace(0,1,n);
    r = 15; % #levellines
    clf; hold on;
    imagesc(t,t,R');
    contour(t,t,R',linspace(0,1,r), 'k');
    colormap(parula(r-1));
    caxis([0 1]);
    plot((X(1:2:end)-1)/n, 'b.', 'MarkerSize', 25);
    plot((X(2:2:end)-1)/n, 'r.', 'MarkerSize', 25);
    axis image; axis off;
    mysaveas('anim2d',it);
    
    % display in 3D
    clf;  hold on;
    surf(t,t,R);
    plot3(imag(X(1:2:end)-1)/n, real(X(1:2:end)-1)/n, X(1:2:end)*0, 'b.', 'MarkerSize', 25);
    plot3(imag(X(2:2:end)-1)/n, real(X(2:2:end)-1)/n, X(2:2:end)*0+1, 'r.', 'MarkerSize', 25);
    view(3);
    shading interp;
    camlight; axis off;
    drawnow;
    mysaveas('anim3d',it);


    
    % advance
    x = x + v;
    % reflexion on boundary
    I = find( real(x)<0 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<0 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
end

