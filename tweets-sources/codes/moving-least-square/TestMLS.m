%%
% Test for MLS

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

p = 10; % #points
d = 0; % degree
% value to interpolate
xi = linspace(0,1,p)';
rand('state', 12);
yi = rand(p,1);

% min sum_i w_i^2 | P(x_i) - y_i  |^2
%    | w_i* |

a = 1;

n = 300;
x = linspace(-.2,1.2,n);
y = MLS(xi,yi,x,d, a);

clf; hold on;
plot(x,y, 'b');
plot(xi,yi, 'r.', 'MarkerSize', 15);
axis([0 1 -.2 1.2]);
box on;
set(gca, 'XTick', [], 'XTick', []);



% animate points in a square with bouncing
rand('state', 123); randn('state', 123);
k = 50; z = rand(p,1)+1i*rand(p,1);
eta = .02; v = randn(p,1) + 1i*randn(p,1); v = eta*v./abs(v);
q = 100;
for it=1:q
    %
    xi = real(z); yi = imag(z);
    y = MLS(xi,yi,x,d, a);
    clf; hold on;
    plot(x,y, 'b', 'LineWidth', 2);
    plot(xi,yi, 'r.', 'MarkerSize', 20);
    axis equal; axis([-.2 1.2 -.2 1.2]);  
    box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(it);
    % DO HERE STUFF
    z = z + v;
    I = find( real(z)<0 | real(z)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(z)<0 | imag(z)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
end


