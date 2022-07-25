addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);



% animate points in a square with bouncing
rand('state', 123); randn('state', 123);
k = 100; x = .1 * ( rand(k,1)+1i*rand(k,1) );
eta = .02; v = randn(k,1) + 1i*randn(k,1); v = eta*v./abs(v);
q = 120;
for it=1:q
    x = x + v;
    I = find( real(x)<0 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<0 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
    
    X = [real(x-mean(x)), imag(x-mean(x))];
    C = (X'*X)/k;
    [U,D] = eig(C); a = sqrt(D(1,1)); b = sqrt(D(2,2));
    u = U(1,1) + 1i*U(2,1);
    w = U(1,2) + 1i*U(2,2);
    t = linspace(0,2*pi,100);
    rho = 1.5;
    
    
    clf;  hold on;
    plot(x, 'k.', 'MarkerSize', 10);
    plot(mean(x), 'r.', 'MarkerSize', 30);
    %
    z = mean(x) +  a*rho*cos(t)*u + b*rho*sin(t)*w;
    h = fill(real(z), imag(z) , 'b', 'LineWidth', 2, 'EdgeColor', 'b');
    h.FaceAlpha = 0.3;
    %
    axis equal; 
    m = .4;
    axis([-m 1+m -m 1+m]); 
    box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(it);
end


