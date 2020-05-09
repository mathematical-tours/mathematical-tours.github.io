% harmonic functions 2D

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep  name '-' znum2str(it,3) '.png']);


n = 250;

rand('state', 123);
randn('state', 123);
% points indicating constraints
k = 6;
x = rand(k,1)+1i*rand(k,1);
a = (-1).^(1:k); a = a(:);
% speed 
eta = .02;
v = randn(k,1) + 1i*randn(k,1);
v = eta*v./abs(v);

proj = @(u)clamp(real(u),1,n) + 1i * clamp(imag(u),1,n);
swap = @(u)imag(u) + 1i*real(u);

phi = @(x)exp(-x);
phi = @(x)x.^2 .* log(x+1e-9);
phi = @(x)exp(-(x/.2).^2);
phi = @(x)x;
%
Phi = @(x,y)phi(abs(x-transpose(y)));

thinplate = 1;

% display grid
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
Z = X + 1i*Y;


q = 80;

for it=1:q
    if thinplate==0
        % solve for weights
        w = inv( Phi(x,x) ) * a;
        % interpolate
        f = reshape( Phi(Z(:),x) * w, [n n]);
    else
        % matlab solution
        p = 1; % interpolation
        p = .8;
        p = 1; % least square
        st = tpaps([real(x) imag(x)]',a', p);
        f = reshape( fnval(st, [real(Z(:)) imag(Z(:))]'), [n n]);
    end
        
    % display quantized colormap
    R = (f+1)/2;
    t = linspace(0,1,n);
    r = 15; % #levellines
    clf; hold on;
    imagesc(t,t,R');
    contour(t,t,R',linspace(-.1,1.1,r), 'k');
    colormap(parula(r-1));
    caxis([-.1 1.1]);
    plot(x(1:2:end), 'b.', 'MarkerSize', 25);
    plot(x(2:2:end), 'r.', 'MarkerSize', 25);
    axis image; axis off;
    drawnow;
    mysaveas('anim2d',it);
    
    % display in 3D
    if 1
    clf;  hold on;
    surf(t,t,R);
    plot3(imag(x(1:2:end)), real(x(1:2:end)), x(1:2:end)*0, 'b.', 'MarkerSize', 25);
    plot3(imag(x(2:2:end)), real(x(2:2:end)), x(2:2:end)*0+1, 'r.', 'MarkerSize', 25);
    view(3);
    shading interp;
    axis([0 1 0 1 -.2 1.2])
    camlight; axis off;
    % drawnow;
    mysaveas('anim3d',it);
    end
    
    % advance
    x = x + v;
    % reflexion on boundary
    I = find( real(x)<0 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<0 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
end

