function draw_linkage(x,y)

hold on;
plot(x, '-', 'LineWidth', 2);
plot(x(2:end), 'r.', 'MarkerSize', 20);
plot(real(x(1)), imag(x(1)), 'b.', 'MarkerSize', 25);
plot(y, 'k.', 'MarkerSize', 25);
axis equal; 
axis([-1 1 -1 1]);
axis off;

end