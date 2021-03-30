% iterates linear/geometric/medians
 
% create save repertory
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

X = randn(3,1) + 1i*randn(3,1);
GeomMeanCpx = @(X)prod(abs(X)).^(1/length(X)) .* exp(1i*mean(angle(X)));

GeomMeanCpx = @(X)prod(abs(X)).^(1/length(X)) .* exp(1i*angle(mean(X./abs(X))));


% animate points in a square with bouncing
rand('state', 123); randn('state', 123);
k = 3; x = 2*( rand(k,1)+1i*rand(k,1) ) - 1 - 1i;
eta = .05; v = randn(k,1) + 1i*randn(k,1); v = eta*v./abs(v);
q = 200;
for it=1:q
    % compute arith-geom
    z = x;
    for j=1:50
        z = [mean(z), GeomMeanCpx(z), FermatPoint(z)];
    end
    % display
    clf; hold on;
    plot(x, 'k.', 'MarkerSize', 30);
    plot(mean(x), 'r.', 'MarkerSize', 30);
    plot(GeomMeanCpx(x), 'b.', 'MarkerSize', 30);
    plot(FermatPoint(x), '.', 'color', [0 .7 0], 'MarkerSize', 30);
    plot(z(1), '.', 'color', [1 .8 0], 'MarkerSize', 30);
    plot([-1 -1 1 1 -1], [-1 1 1 -1 -1], 'k');
    plot([-1 +1], [0 0], 'k:');
    plot([0 0], [-1 +1], 'k:');
    axis equal; axis off; axis([-1.05 1.05 -1.05 1.05]); 
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    % update
    x = x + v;
    I = find( real(x)<-1 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<-1 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
    mysaveas(it);
end
