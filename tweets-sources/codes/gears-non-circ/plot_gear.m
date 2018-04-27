function plot_gear(x, c, color, reverse, t, draw_inside)

% plot_gear - display a gear polar function
%
%    plot_gear(x, c, color, reverse);
%
%   c is the center of the gear, that should be [0 0] for the primal, 
%       and [L 0] for the dual.
%   color is a matlab color for the display.
%   Put reverse=0 for the primal gear, and reverse=1 for the dual geat.
%
%   Copyright (c) 2010 Gabriel Peyre

if nargin<2
    c = [0 0];
end
if nargin<3
    color = 'b';
end
if nargin<4
    reverse = 0;
end
if nargin<6
    draw_inside = 1;
end

c = c(1) + 1i * c(2);

% number of radial line
R = 40;

n = length(x);

lwi = 1;
lw = 2;


if nargin<5 || isempty(t)
    t = 2*pi* (0:n-1)'/n;
end

u = x([1:end 1]) .* exp(1i*t([1:end 1]));
if reverse
    u = -real(u)+1i*imag(u);
end
u = u + c;

hold on;
h = plot( u, color );
set(h, 'LineWidth', lw);
if draw_inside
I = round(linspace(1,n+1,R+1)); I = I(1:end-1);
for i=1:R
    h = plot([c u(I(i))+1e-9i], color);
    set(h, 'LineWidth', lwi);
end
end

h = plot( real(c), imag(c), 'k.' );
set(h, 'MarkerSize', 20);

axis equal;
axis off;