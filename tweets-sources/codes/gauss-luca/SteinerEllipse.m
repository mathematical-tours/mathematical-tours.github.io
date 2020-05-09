%%
% Display for a triangle steiner ellipse.

% create save repertory
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


t = linspace(0,1,512*2);
[Y,X] = meshgrid(t,t);
Z = X + 1i*Y;


% animate points in a square with bouncing
rand('state', 123); randn('state', 123);
k = 3; x = rand(k,1)+1i*rand(k,1);
eta = .02; v = randn(k,1) + 1i*randn(k,1); v = eta*v./abs(v);
q = 120;
for it=1:q
    P = poly(x);
    P1 = polyder(P); x1 = roots(P1);
    P2 = polyder(P1); x2 = roots(P2);
    % mid points
    M = ( x + x([2 3 1]) ) /2;
    U = abs(Z-x1(1)) + abs(Z-x1(2)) <= abs(M(1)-x1(1)) + abs(M(1)-x1(2));
    % draw
    clf; hold on;
    imagesc(t,t, 1-.5*transpose(U));
    colormap gray(256); caxis([0 1]);
    plot(x([1:end 1]), 'r.-', 'LineWidth', 2, 'MarkerSize', 20);
    plot(x1, 'b.', 'MarkerSize', 25);
    % plot(M, '.', 'color', [1 1 1]*.5 + [0 0 1]*.5, 'MarkerSize', 20);
    plot(x2, 'g.', 'MarkerSize', 25);
    axis equal; axis([0 1 0 1]);
    set(gca, 'XTick', [], 'YTick', []); box on;
    drawnow;
    mysaveas(it);
    %
    x = x + v;
    I = find( real(x)<0 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<0 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
end