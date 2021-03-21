

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

addpath('mexEMD/');

a = ones(k,1)/k;
b = ones(k,1)/k;

a = rand(k,1); 
b = rand(k,1); 
a = (1:k)'/k;
b = (1:k)'/k;
a = a/sum(a); b = b/sum(b);


% animate points in a square with bouncing
rand('state', 123); randn('state', 123);
k = 16; x = rand(k,2)+1i*rand(k,2);
x(:,1) = x(:,1)/4;
x(:,2) = x(:,2)/4 + (.75+.75i);
eta = .012; v = randn(k,2) + 1i*randn(k,2); v = eta*v./abs(v);
q = 150;
for it=1:q
    C = abs( x(:,1) - transpose(x(:,2)) );
    [cost,gamma] = mexEMD(a,b,C);
    [I,J,gammaij] = find(gamma);
    %%
    clf; hold on;
    for i=1:length(I)
        plot( [x(I(i),1) x(J(i),2)], 'k' );
    end    
    scatter( real(x(:,1)), imag(x(:,1)), a*k*80, 'b', 'filled' );
    scatter( real(x(:,2)), imag(x(:,2)), b*k*80, 'r', 'filled' );
%    plot(x(:,1), 'r.', 'MarkerSize', 20);
%    plot(x(:,2), 'b.', 'MarkerSize', 20);
    axis equal; axis([0 1 0 1]); box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    % DO HERE STUFF
    x = x + v;
    I = find( real(x)<0 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<0 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
    mysaveas(it);
end