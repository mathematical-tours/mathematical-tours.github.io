%%
% test for finite difference laplacian

N = 128;

x = -N/2:N/2;
clf; hold on;
plot( x/N, 4*sin(pi*x/N).^2, 'r', 'LineWidth', 2 );
plot( x/N, (2*pi*x/N).^2, 'b', 'LineWidth', 2 );
axis tight; box on;