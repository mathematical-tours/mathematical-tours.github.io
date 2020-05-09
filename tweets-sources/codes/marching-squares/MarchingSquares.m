% create save repertory
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

method = 'mc';
method = 'nearest';
method = 'linear';
method = 'spline';

n = 12; 
X = rand(n);


s = .12;
phi = @(d)exp(-d/(2*s^2));
T = .3;
k = 4;
w = ones(k,1); % size

% animate points in a square with bouncing
rand('state', 123); randn('state', 123);
x = rand(k,1)+1i*rand(k,1);
eta = .02; v = randn(k,1) + 1i*randn(k,1); v = eta*v./abs(v);
q = 120;
for it=1:q
    F = GenBall2D([real(x),imag(x)],w,n, phi, T);
    F = rescale(F);
    %
    ls = [1 1]*.5;
    ls = linspace(0,1,6); % ls([1 end]) = [];
    clf; hold on;
    switch method
        case 'mc'
            contourf(1:n,1:n,1-F, ls, 'k', 'LineWidth', 2);
        otherwise
            ups = 30; % upsampling
            s = 1:1/ups:n;
            [Y,X] = meshgrid(1:n,1:n);
            [Ys,Xs] = meshgrid(s,s);
            Fs = interp2(Y,X,F,Ys,Xs, method);
            contourf(s,s,max(1-Fs,0), ls, 'k', 'LineWidth', 2);            
    end
    colormap gray(256);
    plot([1:n;1:n],[ones(1,n);ones(1,n)*n], 'k');
    plot([ones(1,n);ones(1,n)*n], [1:n;1:n],'k');
    axis equal;
    axis([.8 n+.2 0.8 n+.2]);
    caxis([-.3 1]);
    axis off;
    drawnow;
    mysaveas(it);
    %
    x = x + v;
    I = find( real(x)<0 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<0 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
end
