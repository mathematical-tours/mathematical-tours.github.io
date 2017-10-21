%%
% Heat equation on polynomials.

rep = 'results/wave/';
[~,~] = mkdir(rep);
if not(exist('test'))
    test=1;
end

% second order differential on polynomials.
D2 = @(N)diag((1:N-2).*(2:N-1),+2); 

% initial data

%%% from roots

% click and play
clf; hold on;
v = [];
while true
    axis([-1 1 -1 1]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);   
    if button==3
        break;
    end
    v(end+1) = a+1i*b;
end

% click and play
clf; hold on;
w = [];
while true
    axis([-1 1 -1 1]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);   
    if button==3
        break;
    end
    w(end+1) = a+1i*b;
end

tmax = 2;

x= sym('x');
P0 = double( coeffs(expand(prod(x-v)))' );
Q0 = double( coeffs(expand(prod(x-w)))' );


N = max(length(P0),length(Q0));
P0(end+1:N) = 0;
Q0(end+1:N) = 0;

myplot = @(r,m,col)plot(real(r), imag(r), '.', 'MarkerSize', m, 'Color', col);

% use simple explicit euler scheme.
niter = 500; 
tau = tmax/niter;
%
Po = P0;
P  = P0 + tau*Q0;   
clf;
hold on;
ms = 15; B = 1.5;
for i=1:niter
    t = (i-1)/(niter-1);
    Psvg = P;
    P = 2*P - Po + tau^2*D2(N)*P;
    Po = Psvg; 
    %
    myplot( roots(P(end:-1:1)), ms, [t 0 1-t]);
    %
    axis([-B B -B B]); axis equal; box on;
    if mod(i,20)==1
        drawnow;
    end
end
myplot( roots(P0(end:-1:1)), 2*ms, 'k');

axis equal; axis([-1 1 -2 2]); box on;
set(gca, 'FontSize', 15)
saveas(gcf, [rep 'text-' num2str(test) '.eps'], 'epsc');
test = test+1;