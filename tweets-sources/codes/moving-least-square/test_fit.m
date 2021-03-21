%%
% Test for MLS

p = 10;
d = 3; % degree
% value to interpolate
xi = linspace(0,1,p)';
yi = rand(p,1);

% min sum_i w_i^2 | P(x_i) - y_i  |^2
%    | w_i* |

x0 = .3;
s = .06;
w = exp( -(xi-x0).^2 / (2*s^2) );

P = polyfitweighted(xi,yi,d,w);
n = 200;
x = linspace(0,1,n);
y = polyval(P,x);

clf; hold on;
plot(x,y, 'b');
plot(xi,w, 'k--');
plot(xi,yi, 'r.', 'MarkerSize', 15);
axis([0 1 -.2 1.2]);
