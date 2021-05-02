addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep  'anim-' znum2str(it,3) '.png']);



Moments = @(x,y)[x.^3, x.^2.*y,  x.*y.^3, y.^3, x.^2, x.*y, y.^2, x, y  ];
k = 9; 


Moments = @(x,y)[x.^2, x.*y, y.^2, x, y  ];
k = 5; 

Moments = @(x,y)[x.^4, x.^3.*y, x.^2.*y.^2, x.*y.^3, y.^4, ...      
    x.^3, x.^2.*y,  x.*y.^3, y.^3, x.^2, x.*y, y.^2, x, y  ];
k = 14; 


n = 256;
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);


% animate points in a square with bouncing
rand('state', 123); randn('state', 123);

z = rand(k,1)+1i*rand(k,1);
eta = .01*.6; v = randn(k,1) + 1i*randn(k,1); v = eta*v./abs(v);
q = 500;
for it=1:q
    x = real(z); y = imag(z);
    u = pinv( Moments(x,y) ) * ones(k,1);
    F = reshape( Moments(X(:), Y(:)) * u, [n n]);
    %
    clf; hold on;
    contour(t,t,F',[1 1], 'k', 'LineWidth', 3);
    plot(x, y, 'r.', 'MarkerSize', 35);
    axis equal; axis([0 1 0 1]); box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    % DO HERE STUFF
    z = z + v;
    I = find( real(z)<0 | real(z)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(z)<0 | imag(z)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
    mysaveas(it);
end