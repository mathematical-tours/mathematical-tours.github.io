%%
% Display Lagrange interpolation.

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


addpath('../toolbox/');
rep = '../results/lagrange-hermite/';
[~,~] = mkdir(rep);

N = 1024; % discretization
n = 4; % # points
x = linspace(0,1,N);
Q = cumsum(rand(n,1)+.3); Q = rescale(Q,.1,.9);

P = Lagrange(Q,x);

C = distinguishable_colors(n);
clf; hold on;
for i=1:n
    plot(x,P(:,i), 'Color', C(i,:), 'LineWidth', 2);
    plot(Q(i), 1, '.', 'Color', C(i,:), 'MarkerSize', 25);
end
plot(Q,Q*0, '.k', 'MarkerSize', 25);
box on; set(gca, 'XTick', [], 'YTick', []);
axis([0,1,-.5,1.5]);
SetAR(1/2);
saveas(gcf, [rep 'lagrange.eps'], 'epsc');

%%

R = Hermite(Q,x);

clf; hold on;
for i=1:n
    plot(x,R(:,i), 'Color', C(i,:), 'LineWidth', 2);
    plot(Q(i), 1, '.', 'Color', C(i,:), 'MarkerSize', 25);
end
plot(Q,Q*0, '.k', 'MarkerSize', 25);
box on; set(gca, 'XTick', [], 'YTick', []);
axis([0,1,-.5,1.5]);
SetAR(1/2);
saveas(gcf, [rep 'hermite-1.eps'], 'epsc');

clf; hold on;
for i=1:n
    plot(x,R(:,i+n), 'Color', C(i,:), 'LineWidth', 2);
    plot(Q(i), 0, '.', 'Color', C(i,:), 'MarkerSize', 25);
end
box on; set(gca, 'XTick', [], 'YTick', []);
axis([0,1,-.2,.2]);
SetAR(1/2);
saveas(gcf, [rep 'hermite-2.eps'], 'epsc');
