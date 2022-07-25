

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

addpath('mexEMD/');

ka = 8;
kb = 14;

a = ones(ka,1)/ka;
b = ones(kb,1)/kb;

if 0
a = rand(k,1); 
b = rand(k,1); 
a = (1:k)'/k;
b = (1:k)'/k;
end

a = a/sum(a); b = b/sum(b);


% animate points in a square with bouncing
rand('state', 123); randn('state', 123);
xa = rand(ka,1)+1i*rand(ka,1);
xb = rand(kb,1)+1i*rand(kb,1);
xa = xa/4;
xb = xb/4 + (.75+.75i);
eta = .012; 
va = randn(ka,1) + 1i*randn(ka,1); va = eta*va./abs(va);
vb = randn(kb,1) + 1i*randn(kb,1); vb = eta*vb./abs(vb);
q = 300;
for it=1:q
    C = abs( xa - transpose(xb) );
    [cost,gamma] = mexEMD(a,b,C);
    [I,J,gammaij] = find(gamma);
    %%
    clf; hold on;
    for i=1:length(I)
        plot( [xa(I(i)) xb(J(i))], 'k' );
    end    
    scatter( real(xa), imag(xa), a*ka*80, 'b', 'filled' );
    scatter( real(xb), imag(xb), b*kb*80, 'r', 'filled' );
%    plot(x(:,1), 'r.', 'MarkerSize', 20);
%    plot(x(:,2), 'b.', 'MarkerSize', 20);
    s = .02;
    axis equal; axis([-s 1+s -s 1+s]); box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    % DO HERE STUFF
    xa = xa + va;
    xb = xb + vb;
    %
    I = find( real(xa)<0 | real(xa)>1 );
    va(I) = -real(va(I)) + 1i*imag(va(I));
    I = find(  imag(xa)<0 | imag(xa)>1 );
    va(I) = real(va(I)) - 1i*imag(va(I));
    %
    I = find( real(xb)<0 | real(xb)>1 );
    vb(I) = -real(vb(I)) + 1i*imag(vb(I));
    I = find(  imag(xb)<0 | imag(xb)>1 );
    vb(I) = real(vb(I)) - 1i*imag(vb(I));
    mysaveas(it);
end