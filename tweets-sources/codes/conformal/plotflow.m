function plotflow(H1,H2,G1,G2,U)


lw = 2;

clf; hold on;
plot(H1, 'b', 'LineWidth', lw);
plot(H2, 'b', 'LineWidth', lw);
plot(G1, 'r', 'LineWidth', lw);
plot(G2, 'r', 'LineWidth', lw);
fill(real(U),imag(U), [1 1 1]*.8, 'EdgeColor', 'b', 'LineWidth', lw);
axis tight; axis equal; box on; axis off;
% set(gca, 'FontSize', fs);
drawnow;

end