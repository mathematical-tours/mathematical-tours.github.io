% convex hull and delaunay 

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


rand('state', 123);
randn('state', 123);

% points
k = 30; x = rand(k,1)+1i*rand(k,1);
% speed 
eta = .01;
v = randn(k,1) + 1i*randn(k,1); v = eta*v./abs(v);

% to display the paraboloid
c = .5+.5i; d = 1.5; 
xmax = .1;
zmin = -.15*d;

t = linspace(-xmax,1+xmax,500);
[Y,X] = meshgrid(t,t); Z = X + 1i*Y;
Z(abs(Z-c)>.5+xmax) = NaN; 
deltaZ = .02;

q = 120;

for it=1:q
    s = (it-1)/(q-1);
    T = delaunay(real(x), imag(x)); 
   	clf; hold on;  
    if 0
        % planar
        triplot( T, real(x), imag(x), 'r.-', 'LineWidth', 2, 'MarkerSize', 25 );
        axis equal; axis([0 1 0 1]);
    else
        % plot projected on parabola
        t = linspace(0,1,50);
        for i=1:size(T,1)
            lw = 2; 
            u = x(T(i,1))*t + x(T(i,2))*(1-t);
            plot3(real(u), imag(u), d * abs(u-c).^2-deltaZ, 'k', 'LineWidth', lw);
            u = x(T(i,2))*t + x(T(i,3))*(1-t);
            plot3(real(u), imag(u), d * abs(u-c).^2-deltaZ, 'k', 'LineWidth', lw);
            u = x(T(i,3))*t + x(T(i,1))*(1-t);
            plot3(real(u), imag(u), d * abs(u-c).^2-deltaZ, 'k', 'LineWidth', lw);
        end
        ms = 5;
        plot3(real(x), imag(x), d * abs(x-c).^2-deltaZ, 'k.', 'MarkerSize', ms); 
        surf(X,Y, d * abs(Z-c).^2-deltaZ, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
        alpha(.9);
        % plot on the ground
        trimesh(T,real(x), imag(x), zmin + real(x)*0, 'EdgeColor', 'k', 'LineWidth', 1, 'FaceColor', 'none');
        %
        view(3);% view(-60,7);
        axis equal; 
        axis([-xmax, 1+xmax, -xmax, 1+xmax, zmin, .5*d]);
        camlight;
    end
    axis on; box on;
	set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
	drawnow;
    mysaveas(it);   
    % advance
    x = x + v;
    % reflexion on boundary
    I = find( real(x)<0 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<0 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
end
