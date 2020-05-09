% harmonic functions 2D

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep  name '-' znum2str(it,3) '.png']);

% for display purpose
n = 250;

rand('state', 123);
randn('state', 123);

% roots 
k = 15;
x = rand(k,1)+1i*rand(k,1);

% speed 
eta = .02;
v = randn(k,1) + 1i*randn(k,1);
v = eta*v./abs(v);

q = 120;

for it=1:q

    z = x;
    clf;  hold on;
    dmax = min(20,k-2); 
    % dmax = 2;
    for d=1:dmax
        s = (d-1)/(dmax-1);
        ch = convhull(real(z),imag(z));
        plot(z, '.', 'color', [s 0 1-s], 'MarkerSize', 25);
        plot(z(ch([1:end 1])), '-', 'color', [s 0 1-s], 'LineWidth', 2);
        % differentiate polynomials
        P = poly(z);
        % P = (1:(length(P)-1)).*P(2:end);
        P = ((length(P)-1):-1:1).*P(1:end-1);
        z = roots(P);
    end    
    axis equal; box on;
    e = .02;
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

return;

n = 5; % degree

x = sort(randn(5,1));
P = poly(x);
D = randn(size(P));
epsilon = 1E-3;
sort(roots(P+epsilon*D))