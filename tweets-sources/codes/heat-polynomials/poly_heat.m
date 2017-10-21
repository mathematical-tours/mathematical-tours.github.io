%%
% Heat equation on polynomials.

rep = 'results/';
[~,~] = mkdir(rep);
if not(exist('test'))
    test=1;
end

% second order differential on polynomials.
D2 = @(N)diag((1:N-2).*(2:N-1),+2); 

% initial data



P0 = [1 0 -1]'; % X^2-1;

P0 = [1 0 0 0 0 -1]'; % X^5-1;


P0 = [1 0 -1]'; % X^2-1;
P0 = [1 0 0 -1]'; % X^3-1;
P0 = [1 0 0 0 -1]'; % X^3-1;

%%% from roots

% random
v = 2*rand(6,1)-1;
tmax = .5;
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
tmax = .5;

x= sym('x');
P0 = double( coeffs(expand(prod(x-v)))' );


N = length(P0);

myplot = @(r,m,col)plot(real(r), imag(r), '.', 'MarkerSize', m, 'Color', col);

% use simple explicit euler scheme.
niter = 500; 
tau = tmax/niter;
%
P = P0;
clf;
hold on;
ms = 15; B = 1.5;
for i=1:niter
    t = (i-1)/(niter-1);
    P = P + tau*D2(N)*P;
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