% convex hull in 2D

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep  name '-' znum2str(it,3) '.png']);


rand('state', 123);
randn('state', 123);

% roots 
k = 50;
x = rand(k,1)+1i*rand(k,1);
% speed 
eta = .02;
v = randn(k,1) + 1i*randn(k,1);
v = eta*v./abs(v);

q = 120;

t = linspace(-100,100,1000);

for it=1:q
    s = (it-1)/(q-1);
    %
    ch = convhull(real(x),imag(x));
    %
    clf;  hold on;
    pgon = polyshape(real(x(ch)),imag(x(ch)));
    plot(pgon,'FaceColor', [s 0 1-s],'FaceAlpha',0.1);
    for i=1:length(ch)
        j = mod(i, length(ch))+1;
        plot( t*x(ch(i)) + (1-t)*x(ch(j)), 'k' );
    end
    plot(x, '.', 'color', [s 0 1-s], 'MarkerSize', 10);
    plot(x(ch), '.', 'color', [s 0 1-s], 'MarkerSize', 30);
    plot(x(ch([1:end 1])), '-', 'color', [s 0 1-s], 'LineWidth', 2);
    axis equal; box on;
    e = .15;
    axis([-e 1+e -e 1+e]);
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas('anim',it);
   
    % advance
    x = x + v;
    % reflexion on boundary
    I = find( real(x)<0 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<0 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
end
